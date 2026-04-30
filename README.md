# ScaffoldTuner

A Python toolkit for proposing structurally modified drug-like molecules based on Murcko scaffolds.  
Given a parent molecule (SMILES), ScaffoldTuner systematically replaces side chains to increase or decrease specific molecular features such as hydrogen bond donors/acceptors, aromatic rings, rotatable bonds, and heteroatom counts.

---

## Features

- Extracts Murcko scaffold from a parent molecule
- Decomposes R-groups using RDKit's R-group decomposition
- Proposes candidate molecules by replacing one substituent at a time
- Supports 10 modification modes (see below)
- Filters candidates by whether the target feature changes in the intended direction
- Optional `strict` mode: ensures only the target feature changes

---

## Requirements

- Python 3.8+
- [RDKit](https://www.rdkit.org/)

```bash
conda install -c conda-forge rdkit
```

---

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/scaffold-tuner.git
cd scaffold-tuner
```

---

## Usage

```python
from scaffold_tuner.scaffold_intervention import propose_structures, print_proposals

parent_smiles = "c1ccc(NC(=O)c2ccccc2)cc1"  # Example: benzanilide

results = propose_structures(
    parent_smiles=parent_smiles,
    mode="add_hbd",        # Increase hydrogen bond donors
    max_candidates=5,
    random_seed=42,
)

print_proposals(results)
```

### Example Output

```
=== Candidate 1 ===
mode: add_hbd,  site: R1,  fragment: [*:1]N
parent SMILES   : O=C(Nc1ccccc1)c1ccccc1
generated SMILES: Nc1ccc(NC(=O)c2ccccc2)cc1
features ['NumHAcceptors', 'NumHDonors', 'NumAromaticRings', 'NumRotatableBonds', 'NumHeteroatoms']
  parent   : (1, 1, 2, 2, 1)
  generated: (2, 2, 2, 2, 2)
```

---

## Modification Modes

| Mode | Description |
|---|---|
| `add_hbd` | Increase hydrogen bond donors |
| `remove_hbd` | Decrease hydrogen bond donors |
| `add_hba` | Increase hydrogen bond acceptors |
| `remove_hba` | Decrease hydrogen bond acceptors |
| `add_ar` | Increase number of aromatic rings |
| `remove_ar` | Decrease number of aromatic rings |
| `add_rb` | Increase number of rotatable bonds |
| `remove_rb` | Decrease number of rotatable bonds |
| `add_ha` | Increase number of heteroatoms |
| `remove_ha` | Decrease number of heteroatoms |

---

## Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `parent_smiles` | str | required | SMILES of the parent molecule |
| `mode` | str | `"add_hbd"` | Modification mode (see table above) |
| `max_cuts` | int or None | None | Max number of substitution sites |
| `max_candidates` | int | 10 | Maximum number of candidates to return |
| `shuffle_sites` | bool | False | Randomly shuffle substitution sites |
| `shuffle_fragments` | bool | False | Randomly shuffle fragment candidates |
| `random_seed` | int or None | None | Random seed for reproducibility |
| `strict` | bool | False | If True, only return candidates where no other features change |

---

## Project Structure

```
scaffold-tuner/
├── README.md
├── LICENSE
├── .gitignore
└── scaffold_tuner/
    ├── __init__.py
    └── scaffold_intervention.py
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
