from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import  BondType, rdFMCS
import numpy as np
from tqdm import tqdm
from collections import defaultdict, namedtuple
import itertools
from sklearn.cluster import DBSCAN


def rdkit_MCS_Sim(mol1, mol2):
    res = rdFMCS.FindMCS([mol1, mol2], timeout = 1, ringMatchesRingOnly=True,completeRingsOnly=True)
    return res, (2*res.numAtoms)/(mol1.GetNumAtoms()+mol2.GetNumAtoms())

def murcko_input(smiles_list):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
    murcko_mols = [MurckoScaffold.GetScaffoldForMol(m) for m in mols]
    murcko_smiles_list = list(set([Chem.MolToSmiles(s) for s in murcko_mols]))
    return murcko_smiles_list

def make_csk(smiles):
    mol = Chem.MolFromSmiles(smiles)
    murcko_mol = MurckoScaffold.GetScaffoldForMol(mol)
    csk_mol = MurckoScaffold.MakeScaffoldGeneric(murcko_mol)
    return {mol: csk_mol}, {csk_mol: Chem.MolToSmiles(csk_mol)}

def find_key_by_value(dictionary, value):
    for k, v in dictionary.items():
        if v == value:
            return k

def smiles_list_processing(smiles_list):
    murcko_smiles_list = murcko_input(smiles_list)
    mol_csk_data = [make_csk(smi) for smi in murcko_smiles_list]
    mol_csk_dict = {list(d[0].keys())[0] : list(d[0].values())[0] for d in mol_csk_data}
    csk_smi_dict = {list(d[1].keys())[0] : list(d[1].values())[0] for d in mol_csk_data}
    csk_smiles_list = list(set(smi for smi in csk_smi_dict.values()))
    csk_mols = [find_key_by_value(csk_smi_dict, smi) for smi in csk_smiles_list]
    csk_sim_matrix = np.array([[rdkit_MCS_Sim(m1, m2)[1] for m2 in csk_mols] for m1 in csk_mols])


    dist = 1.0 - csk_sim_matrix
    np.fill_diagonal(dist, 0.0)

    db = DBSCAN(eps=0.1, min_samples=1, metric="precomputed")
    labels = db.fit_predict(dist)

    return (mol_csk_dict, csk_smi_dict), (csk_mols, labels)

def dedup_keep_max(mol_sim_dict):
    """Collapse structurally identical molecules, keeping the highest value.
    Returns a new dict {representative_mol: max_similarity}."""
    best = {}
    for mol, val in mol_sim_dict.items():
        key = None
        for fn in (Chem.MolToSmiles, Chem.MolToSmarts):
            try:
                s = fn(mol)
                if s:
                    key = s
                    break
            except Exception:
                pass
        if key is None:
            key = id(mol)
        if key not in best or val > best[key][1]:
            best[key] = (mol, val)
    return {mol: val for mol, val in best.values()}

def has_broken_ring(scaffold):
    Chem.GetSSSR(scaffold)
    ring_info = scaffold.GetRingInfo()
    for bond in scaffold.GetBonds():
        if bond.HasQuery():
            dq = bond.DescribeQuery()
            requires_ring = ("BondInRing 1 = val" in dq
                             or bond.GetBondType() == Chem.BondType.AROMATIC)
        else:
            requires_ring = bond.IsInRing()
        if requires_ring and not ring_info.NumBondRings(bond.GetIdx()):
            return True
    return False

def mcs_murcko(mol):
    """Murcko-style scaffold that also works on MCS query mols.
    Iteratively deletes terminal (degree <= 1) atoms that aren't in a ring,
    leaving only ring systems and the linkers connecting them."""
    rw = Chem.RWMol(mol)
    while True:
        Chem.GetSSSR(rw)
        ri = rw.GetRingInfo()
        leaves = [a.GetIdx() for a in rw.GetAtoms()
                  if a.GetDegree() <= 1 and not ri.NumAtomRings(a.GetIdx())]
        if not leaves:
            break
        for idx in sorted(leaves, reverse=True):
            rw.RemoveAtom(idx)
    return rw.GetMol()


def build_clusters(cmols, labels):
	clusters = dict()
	for k in tqdm(sorted(set(labels))):
		query_dict = dict()
		idx = [i for i in range(len(cmols)) if labels[i] == k]
		mols_k = [cmols[i] for i in idx]

		if len(mols_k) == 1:
			clusters[k] = {mols_k[0]: 1}
			continue

		n = len(mols_k)
		for i in tqdm(range(n)):
			for j in range(i + 1, n):
				res, mcs_sim = rdkit_MCS_Sim(mols_k[i], mols_k[j])
				if mcs_sim < 0.9:
					continue
				if has_broken_ring(res.queryMol):
					continue
				new_query_mol = mcs_murcko(res.queryMol)
				query_dict[new_query_mol] = mcs_sim
		clusters[k] = dedup_keep_max(query_dict)
	return clusters


def refine_cluster(clusters, k=0, n_iterations=5):
	for iteration in tqdm(range(n_iterations), desc="iterations"):
		used_scaffolds = []

		if len(clusters[k]) == 1:
			break

		prev = dict(clusters[k])
		new_scaffolds = dict()
		for comb in itertools.combinations(clusters[k], 2):

			if any(sc in used_scaffolds for sc in comb):
				continue

			parents_sims = clusters[k][comb[0]], clusters[k][comb[1]]
			res, mcs_sim = rdkit_MCS_Sim(comb[0], comb[1])
			if mcs_sim >= max(*parents_sims):
				if has_broken_ring(res.queryMol):
					continue
				else:
					new_query_mol = mcs_murcko(res.queryMol)
					new_scaffolds[new_query_mol] = mcs_sim
					[used_scaffolds.append(sc) for sc in comb if sc not in used_scaffolds]
			else:
				new_scaffolds[comb[0]] = clusters[k][comb[0]]
				new_scaffolds[comb[1]] = clusters[k][comb[1]]
		clusters[k] = dedup_keep_max(new_scaffolds)
		if clusters[k] == prev:
			break
	return clusters


bond_types = {
	Chem.BondType.SINGLE: '-',
	Chem.BondType.DOUBLE: '=',
	Chem.BondType.TRIPLE: '#',
	Chem.BondType.AROMATIC: ':',
	Chem.BondType.UNSPECIFIED: '~',
}


def get_bonds_from_match(query_mol, atom_match):
	cind_mind_bond_dict = dict()
	for bond in query_mol.GetBonds():
		idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
		cind_mind_bond_dict[tuple(sorted((atom_match[idx1], atom_match[idx2])))] = tuple(sorted((idx1, idx2)))
	return cind_mind_bond_dict


def atoms_bonds_properties_extractor(initial_mol, csk_mol, mcs_mol):
	info = namedtuple('info', ['symbol', 'num', 'ring'])
	atom_props = dict()
	bond_props = dict()
	atom_match = csk_mol.GetSubstructMatch(mcs_mol)
	cind_mind_bond_dict = get_bonds_from_match(mcs_mol, atom_match)

	if not atom_match:
		return False

	cind_mind_atom_dict = {ind: i for i, ind in enumerate(atom_match)}

	for aid in range(initial_mol.GetNumAtoms()):
		if aid in cind_mind_atom_dict:
			a = initial_mol.GetAtomWithIdx(aid)
			atom_props[cind_mind_atom_dict[aid]] = info(
				(str.lower(a.GetSymbol()) if a.GetIsAromatic() else a.GetSymbol()), a.GetAtomicNum(), a.IsInRing())

	for b in initial_mol.GetBonds():
		begin, end = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
		bid = tuple(sorted((begin, end)))
		if bid in cind_mind_bond_dict:
			bond_props[cind_mind_bond_dict[bid]] = b.GetBondType()
	return atom_props, bond_props


def record_to_obs(rec):
	num, sym = rec.num, rec.symbol
	extra = []
	if getattr(rec, 'ring', None) is not None:
		extra.append('R' if rec.ring else '!R')
	return num, sym, tuple(extra)


def merge_atom(observations):
	"""observations = list of records for the SAME index across all sources."""
	groups = defaultdict(set)
	for rec in observations:
		num, sym, extra = record_to_obs(rec)
		groups[(num, extra)].add(sym)

	specs = []
	for (num, extra), syms in groups.items():
		elem = f"#{num}" if len(syms) > 1 else next(iter(syms))
		specs.append("&".join([elem, *extra]))
	return "[" + ",".join(sorted(set(specs))) + "]"


def merge_maps(maps):
	"""maps = list of {idx: record}.  -> {idx: smarts}."""
	per_idx = defaultdict(list)
	for m in maps:
		for idx, rec in m.items():
			per_idx[idx].append(rec)
	return {idx: merge_atom(recs) for idx, recs in sorted(per_idx.items())}


def merge_bonds(maps):
	"""maps = list of {bond_idx: BondType}.  -> {bond_idx: BondType}.
	   Keep the type if every source agrees; otherwise UNSPECIFIED."""
	per_idx = defaultdict(set)
	for m in maps:
		for idx, bt in m.items():
			per_idx[idx].add(bt)
	return {idx: (next(iter(bts)) if len(bts) == 1 else BondType.UNSPECIFIED)
	        for idx, bts in sorted(per_idx.items())}


def mol_to_query_mol(mcs_scaffold, merged_atom_queries, merged_bond_queries):
	new_mcs = Chem.RWMol(mcs_scaffold)

	for idx, query in merged_atom_queries.items():
		qa = Chem.AtomFromSmarts(query)
		new_mcs.ReplaceAtom(idx, qa)

	for idx, query in merged_bond_queries.items():
		nb = Chem.BondFromSmarts(bond_types[query])
		idx = new_mcs.GetBondBetweenAtoms(idx[0], idx[1]).GetIdx()
		new_mcs.ReplaceBond(idx, nb)

	return new_mcs


def scaffold_builder(mol_csk_dict, mcs_scaffold):
	atom_queries = list()
	bond_queries = list()
	for mol, csk in mol_csk_dict.items():
		atom_info, bond_info = atoms_bonds_properties_extractor(mol, csk, mcs_scaffold)
		atom_queries.append(atom_info)
		bond_queries.append(bond_info)
	merged_atom_queries = merge_maps(atom_queries)
	merged_bond_queries = merge_bonds(bond_queries)
	scaffold_mol = mol_to_query_mol(mcs_scaffold, merged_atom_queries, merged_bond_queries)
	if scaffold_mol.GetNumAtoms() > 6:
		new_scaffold = Chem.MolToSmarts(mcs_murcko(scaffold_mol))
	else:
		return False
	return new_scaffold


def generate_scaffolds(smiles_list, n_iterations=15):
	mol_data, cluster_data = smiles_list_processing(smiles_list)
	mol_csk_dict, csk_smi_dict = mol_data
	csk_mols, labels = cluster_data

	clusters = build_clusters(csk_mols, labels)
	for k in clusters.keys():
		clusters = refine_cluster(clusters, k=k, n_iterations=n_iterations)

	scaffolds = dict()
	for k in clusters.keys():
		for mcs, score in clusters[k].items():
			matched_mol_csk_dict = {mol: csk for mol, csk in mol_csk_dict.items() if csk.HasSubstructMatch(mcs)}
			general_scaffold = scaffold_builder(matched_mol_csk_dict, mcs)
			if general_scaffold:
				scaffolds[general_scaffold] = (score, len(matched_mol_csk_dict))
	return scaffolds

def coverage(smiles_list, scaffold_list, use_smarts=False):
    parser = Chem.MolFromSmarts if use_smarts else Chem.MolFromSmiles
    scaffs = [parser(s) for s in scaffold_list]
    scaffs = [s for s in scaffs if s is not None]

    covered, uncovered = [], []
    for smi in smiles_list:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            uncovered.append(smi)
            continue
        if any(m.HasSubstructMatch(sc) for sc in scaffs):
            covered.append(smi)
        else:
            uncovered.append(smi)

    pct = 100.0 * len(covered) / len(smiles_list) if smiles_list else 0.0
    return pct, covered, uncovered