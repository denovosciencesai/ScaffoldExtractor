# General Scaffold Extractor

A lightweight tool designed to extract general scaffolds from a list of chemical structures (SMILES) and assess how well these scaffolds cover your initial molecular dataset.

---

## 🚀 Features

* **Scaffold Generation**: Extract core chemical scaffolds from your compound library.
* **Coverage Analysis**: Calculate what percentage of your original molecules are covered by the generated scaffolds.
* **Fast Setup**: Includes a pre-configured Conda environment setup.

---

## 📦 Installation

To get started, make sure you have [Conda](https://docs.conda.io/en/latest/) or [Miniconda](https://docs.conda.io/en/miniconda.html) installed.

### 1. Create the Environment
This repository ships with an `environment.yml`. From the repository root, run:

```bash
conda env create -f environment.yml
```

### 2. Activate the Environment
Once the installation is complete, activate the environment using:
```bash
conda activate scaffold_extractor
```

---

## 💻 Quick Start & Usage

To use the tool, place your input `.txt` file (containing your SMILES, one per line) and the `GeneralScaffoldExtractor.py` file in the same directory as your script.

You can copy, paste, and run this Python code directly:

```python
from GeneralScaffoldExtractor import generate_scaffolds, coverage

# 1. Load your list of molecules (SMILES) from the .txt file
# Replace 'smiles_list.txt' with the actual name of your text file
with open("smiles_list.txt", "r") as f:
    list_of_smiles = [line.strip() for line in f if line.strip()]

# 2. Extract scaffolds
scaffolds = generate_scaffolds(list_of_smiles)
print(f"Generated {len(scaffolds)} scaffolds.")

# 3. Assess the coverage of the generated scaffolds
# This checks what part of the initial molecules is covered by the scaffolds
pct, covered, uncovered = coverage(list_of_smiles, scaffolds, use_smarts=True)
print(f"Coverage: {pct:.2f}% ({len(covered)} covered, {len(uncovered)} uncovered)")
```

Running the above on the bundled `smiles_list.txt` (138 molecules) produces 19 scaffolds
covering 99.28% of the input set.

> ⚠️ **Note:** Always ensure `GeneralScaffoldExtractor.py` and your `.txt` file are in the same folder as your execution script to avoid path-related errors.

---

## 🛠️ API & Function Reference

The `GeneralScaffoldExtractor` module exposes two primary functions to streamline your workflow:

### 1. `generate_scaffolds(smiles_list, n_iterations=15)`
Extracts core structural scaffolds from a collection of molecules.

| Parameter | Type  | Default | Description |
| :--- |:------| :--- | :--- |
| **`smiles_list`** | `list` | *Required* | A list containing SMILES strings of your initial dataset. |
| **`n_iterations`** | `int` | `15` | Maximum number of refinement passes used to merge scaffolds within a cluster. Refinement stops early once a cluster converges. |

* **Returns**: `dict`
  * Maps each generated scaffold (a **SMARTS** string) to a tuple `(score, n_molecules)`, where `score` is the MCS similarity the scaffold was derived from and `n_molecules` is the number of input structures matched by that scaffold.

---

### 2. `coverage(smiles_list, scaffold_list, use_smarts=False)`
Evaluates how effectively the generated scaffolds represent (cover) your original molecular library by performing substructure matches.

| Parameter | Type   | Default | Description |
| :--- |:-------| :--- | :--- |
| **`smiles_list`** | `list` | *Required* | The original list of SMILES molecules you started with. |
| **`scaffold_list`** | `iterable` | *Required* | The scaffolds to match against. The `dict` returned by `generate_scaffolds` can be passed directly, since iterating it yields the scaffold strings. |
| **`use_smarts`** | `bool` | `False` | Parse `scaffold_list` entries as SMARTS rather than SMILES. Set this to `True` for the output of `generate_scaffolds`, which returns SMARTS. |

* **Returns**: `tuple` — `(pct, covered, uncovered)`
  * **`pct`** (`float`): percentage of the input molecules matched by at least one scaffold.
  * **`covered`** (`list`): the input SMILES that matched at least one scaffold.
  * **`uncovered`** (`list`): the input SMILES that matched none. SMILES that cannot be parsed are also reported here.

---

## 📄 License

Released under the [MIT License](LICENSE).
