# -*- coding: utf-8 -*-
import random
from rdkit import Chem
from rdkit.Chem import Lipinski, Descriptors, rdRGroupDecomposition
from rdkit.Chem.Scaffolds import MurckoScaffold


# ──────────────────────────────────────────────
# 1. Parent molecule → Murcko core with R-group labels
# ──────────────────────────────────────────────
def make_murcko_core_with_rlabels(parent_smiles, max_cuts=None):
    """
    Generate the Murcko scaffold from the parent molecule SMILES and
    return the core with attachment point labels [*:1], [*:2], ... at
    each scaffold-sidechain boundary.

    Returns: (parent_mol, scaffold_mol, core_with_labels)
    """
    parent_mol = Chem.MolFromSmiles(parent_smiles)
    if parent_mol is None:
        raise ValueError("Invalid parent SMILES.")

    scaffold_mol = MurckoScaffold.GetScaffoldForMol(parent_mol)
    if scaffold_mol is None:
        raise ValueError("Failed to generate Murcko scaffold.")

    matches = parent_mol.GetSubstructMatches(scaffold_mol)
    if not matches:
        raise ValueError("Murcko scaffold did not match the parent molecule.")

    # Atom indices in the parent molecule that belong to the scaffold,
    # and a reverse mapping from parent index to scaffold index
    match = matches[0]
    core_atom_idxs = set(match)
    parent_to_core = {p: c for c, p in enumerate(match)}

    # Collect bonds that cross the scaffold-sidechain boundary
    boundary_bonds = [
        b for b in parent_mol.GetBonds()
        if (b.GetBeginAtomIdx() in core_atom_idxs) ^ (b.GetEndAtomIdx() in core_atom_idxs)
    ]
    if max_cuts is not None:
        boundary_bonds = boundary_bonds[:max_cuts]

    # For each boundary bond, add a dummy atom [*:N] to the scaffold-side atom
    rw_core = Chem.RWMol(scaffold_mol)
    for rlabel, bond in enumerate(boundary_bonds, start=1):
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        core_parent_idx = a1 if a1 in core_atom_idxs else a2
        core_idx = parent_to_core[core_parent_idx]

        dummy = Chem.Atom(0)          # dummy atom *
        dummy.SetAtomMapNum(rlabel)
        dummy_idx = rw_core.AddAtom(dummy)
        rw_core.AddBond(core_idx, dummy_idx, bond.GetBondType())

    core_with_labels = rw_core.GetMol()
    Chem.SanitizeMol(core_with_labels)
    return parent_mol, scaffold_mol, core_with_labels


# ──────────────────────────────────────────────
# 2. Utilities
# ──────────────────────────────────────────────
def count_features(mol):
    """Return a tuple of (HBA, HBD, NumAromaticRings, NumRotatableBonds, NumHeteroatoms)."""
    return (
        Lipinski.NumHAcceptors(mol),
        Lipinski.NumHDonors(mol),
        Descriptors.NumAromaticRings(mol),
        Descriptors.NumRotatableBonds(mol),
        Descriptors.NumHeteroatoms(mol),
    )


def get_r_labels(core_mol):
    """Return a sorted list of atom map numbers for all dummy atoms [*:N] in the core."""
    return sorted({
        atom.GetAtomMapNum()
        for atom in core_mol.GetAtoms()
        if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() > 0
    })


# ──────────────────────────────────────────────
# 3. Extract the original R-groups from the parent molecule
# ──────────────────────────────────────────────
def get_original_rgroups(parent_smiles, core_with_labels):
    """
    Perform R-group decomposition on the parent molecule and return
    each substituent as {"R1": SMILES, "R2": SMILES, ...}.
    """
    parent_mol = Chem.MolFromSmiles(parent_smiles)
    if parent_mol is None:
        raise ValueError("Invalid parent SMILES.")

    params = rdRGroupDecomposition.RGroupDecompositionParameters()
    params.labels = rdRGroupDecomposition.RGroupLabels.AtomMapLabels
    params.rgroupLabelling = rdRGroupDecomposition.RGroupLabelling.AtomMap
    params.onlyMatchAtRGroups = True

    rgd = rdRGroupDecomposition.RGroupDecomposition(core_with_labels, params)
    if rgd.Add(parent_mol) < 0:
        raise ValueError("RGroupDecomposition.Add() failed.")
    if not rgd.Process():
        raise RuntimeError("RGroupDecomposition.Process() failed.")

    rows = rgd.GetRGroupsAsRows()
    if not rows:
        raise RuntimeError("No R-group decomposition result was obtained.")

    # Exclude the "Core" key and convert each R-group molecule to SMILES
    return {
        key: Chem.MolToSmiles(mol)
        for key, mol in rows[0].items()
        if key != "Core" and mol is not None
    }


# ──────────────────────────────────────────────
# 4. Fragment attachment
# ──────────────────────────────────────────────
def _get_dummy_info(mol, map_num):
    """Find the dummy atom [*:map_num] and return (dummy_idx, neighbor_idx, bond_type)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() == map_num:
            nbr = atom.GetNeighbors()[0]
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            return atom.GetIdx(), nbr.GetIdx(), bond.GetBondType()
    raise ValueError(f"Dummy atom [*:{map_num}] not found.")


def attach_fragment(scaffold_mol, fragment_mol, map_num):
    """Connect the [*:map_num] attachment points of scaffold and fragment, then remove the dummy atoms."""
    s_dummy, s_nbr, s_bond = _get_dummy_info(scaffold_mol, map_num)
    f_dummy, f_nbr, _      = _get_dummy_info(fragment_mol, map_num)

    offset = scaffold_mol.GetNumAtoms()
    emol = Chem.EditableMol(Chem.CombineMols(scaffold_mol, fragment_mol))
    emol.AddBond(s_nbr, f_nbr + offset, order=s_bond)

    # Remove dummy atoms from the higher index first to avoid index shifting
    for idx in sorted([s_dummy, f_dummy + offset], reverse=True):
        emol.RemoveAtom(idx)

    mol = emol.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def build_molecule(scaffold_smiles, rgroups):
    """
    Attach all R-groups in rgroups (e.g., {"R1": "[*:1]N", "R2": "[*:2]OC"})
    to scaffold_smiles in ascending order of R-group number.
    """
    mol = Chem.MolFromSmiles(scaffold_smiles)
    if mol is None:
        raise ValueError("Invalid scaffold SMILES.")

    # Attach in ascending R-group number order to avoid order dependency
    for site in sorted(rgroups, key=lambda x: int(x[1:])):
        frag = Chem.MolFromSmiles(rgroups[site])
        if frag is None:
            raise ValueError(f"Invalid fragment SMILES: {rgroups[site]}")
        mol = attach_fragment(mol, frag, int(site[1:]))
    return mol


# ──────────────────────────────────────────────
# 5. Fragment pools for increasing/decreasing each descriptor
# ──────────────────────────────────────────────

# All templates use [*:1] as the attachment point placeholder.
# At runtime, [*:1] is replaced with the actual R-group label [*:N].
_FRAGMENT_TEMPLATES = {
    # HBD (hydrogen bond donors)
    "add_hbd": [
        "[*:1]N",           # amino
        "[*:1]NC",          # methylamino
        "[*:1]O",           # hydroxy
        "[*:1]NS(=O)(=O)C", # sulfonamide-like
    ],
    "remove_hbd": [
        "[*:1]C",           # methyl (no HBD)
        "[*:1]CC",          # ethyl (no HBD)
        "[*:1]F",           # fluoro (no HBD)
        "[*:1]Cl",          # chloro (no HBD)
        "[*:1]C(C)C",       # isopropyl (no HBD)
    ],
    # HBA (hydrogen bond acceptors)
    "add_hba": [
        "[*:1]F",           # fluoro
        "[*:1]OC",          # methoxy
        "[*:1]C#N",         # cyano
        "[*:1]C(=O)C",      # keto-like
        "[*:1]N(C)C",       # tertiary amine
    ],
    "remove_hba": [
        "[*:1]C",           # methyl (no HBA)
        "[*:1]CC",          # ethyl (no HBA)
        "[*:1]F",           # fluoro (counted as HBA=0 by RDKit in most cases)
        "[*:1]Cl",          # chloro (no HBA)
        "[*:1]c1ccccc1",    # phenyl (no HBA)
    ],
    # NumAromaticRings
    "add_ar": [
        "[*:1]c1ccccc1",    # phenyl
        "[*:1]c1ccncc1",    # pyridyl
        "[*:1]c1ccco1",     # 2-furyl
        "[*:1]c1cccs1",     # 2-thienyl
    ],
    "remove_ar": [
        "[*:1]C",           # methyl (no aromatic ring)
        "[*:1]CC",          # ethyl (no aromatic ring)
        "[*:1]CCC",         # propyl (no aromatic ring)
        "[*:1]C(C)C",       # isopropyl (no aromatic ring)
        "[*:1]C1CCCCC1",    # cyclohexyl (no aromatic ring)
    ],
    # NumRotatableBonds
    "add_rb": [
        "[*:1]CC",          # ethyl (+1 rotatable bond)
        "[*:1]CCC",         # propyl (+2 rotatable bonds)
        "[*:1]CCCC",        # butyl (+3 rotatable bonds)
        "[*:1]COC",         # methoxymethyl (+2 rotatable bonds)
        "[*:1]CC(=O)C",     # acetonyl (+2 rotatable bonds)
    ],
    "remove_rb": [
        "[*:1]C",           # methyl (no rotatable bonds)
        "[*:1]F",           # fluoro (no rotatable bonds)
        "[*:1]Cl",          # chloro (no rotatable bonds)
        "[*:1]C1CC1",       # cyclopropyl (restricted rotation by ring)
        "[*:1]C1CCCCC1",    # cyclohexyl (restricted rotation by ring)
    ],
    # NumHeteroatoms
    "add_ha": [
        "[*:1]N",           # amino (adds N)
        "[*:1]O",           # hydroxy (adds O)
        "[*:1]OC",          # methoxy (adds O)
        "[*:1]F",           # fluoro (adds F)
        "[*:1]S",           # thiol (adds S)
    ],
    "remove_ha": [
        "[*:1]C",           # methyl (no heteroatoms)
        "[*:1]CC",          # ethyl (no heteroatoms)
        "[*:1]C(C)C",       # isopropyl (no heteroatoms)
        "[*:1]c1ccccc1",    # phenyl (no heteroatoms)
        "[*:1]C1CCCCC1",    # cyclohexyl (no heteroatoms)
    ],
}

# mode → (descriptor index, direction)
# Descriptor index corresponds to the return value of count_features():
#   0=HBA, 1=HBD, 2=NumAromaticRings, 3=NumRotatableBonds, 4=NumHeteroatoms
# Direction: +1 = increase expected, -1 = decrease expected
_MODE_CONFIG = {
    "add_hbd":    (1, +1),
    "remove_hbd": (1, -1),
    "add_hba":    (0, +1),
    "remove_hba": (0, -1),
    "add_ar":     (2, +1),
    "remove_ar":  (2, -1),
    "add_rb":     (3, +1),
    "remove_rb":  (3, -1),
    "add_ha":     (4, +1),
    "remove_ha":  (4, -1),
}

def get_fragment_pool(mode, map_num):
    """Return the fragment templates for the given mode with [*:1] replaced by [*:map_num]."""
    if mode not in _FRAGMENT_TEMPLATES:
        raise ValueError(f"mode must be one of {list(_MODE_CONFIG.keys())}")
    return [t.replace("[*:1]", f"[*:{map_num}]") for t in _FRAGMENT_TEMPLATES[mode]]


# ──────────────────────────────────────────────
# 6. Main function: propose candidate structures
# ──────────────────────────────────────────────
def propose_structures(
    parent_smiles,
    mode="add_hbd",
    max_cuts=None,
    max_candidates=10,
    shuffle_sites=False,
    shuffle_fragments=False,
    random_seed=None,
    strict=False,
):
    """
    Propose candidate structures in which one descriptor of the parent molecule
    has been increased or decreased. Unmodified R-group positions retain the
    original substituents from the parent molecule.

    Parameters
    ----------
    parent_smiles     : SMILES of the parent molecule
    mode              : Transformation mode. One of the following:
                          "add_hbd"    : increase hydrogen bond donor count
                          "remove_hbd" : decrease hydrogen bond donor count
                          "add_hba"    : increase hydrogen bond acceptor count
                          "remove_hba" : decrease hydrogen bond acceptor count
                          "add_ar"     : increase aromatic ring count
                          "remove_ar"  : decrease aromatic ring count
                          "add_rb"     : increase rotatable bond count
                          "remove_rb"  : decrease rotatable bond count
                          "add_ha"     : increase heteroatom count
                          "remove_ha"  : decrease heteroatom count
    max_cuts          : Maximum number of attachment points to label on the scaffold
                        (None = use all boundary bonds)
    max_candidates    : Maximum number of candidates to return
    shuffle_sites     : If True, explore attachment sites in random order
    shuffle_fragments : If True, explore fragments in random order
    random_seed       : Random seed for reproducibility
    strict            : If True, return only candidates in which exclusively the
                        target descriptor has changed (all other descriptors unchanged)
    """
    if mode not in _MODE_CONFIG:
        raise ValueError(f"mode must be one of {list(_MODE_CONFIG.keys())}")

    if random_seed is not None:
        random.seed(random_seed)

    parent_mol, _, core = make_murcko_core_with_rlabels(parent_smiles, max_cuts=max_cuts)
    core_smiles = Chem.MolToSmiles(core)
    labels = get_r_labels(core)
    if not labels:
        return []

    parent_canon = Chem.MolToSmiles(parent_mol)
    parent_feats = count_features(parent_mol)
    feat_idx, direction = _MODE_CONFIG[mode]  # descriptor index and direction of change

    original_rgroups = get_original_rgroups(parent_smiles, core)

    if shuffle_sites:
        labels = random.sample(labels, len(labels))

    results, seen = [], set()

    for label in labels:
        site = f"R{label}"
        pool = get_fragment_pool(mode, label)
        if shuffle_fragments:
            pool = random.sample(pool, len(pool))

        for frag_smiles in pool:
            try:
                # Replace only one R-group and build the new molecule
                rgroups = {**original_rgroups, site: frag_smiles}
                gen_mol = build_molecule(core_smiles, rgroups)
                gen_canon = Chem.MolToSmiles(gen_mol)

                # Skip if identical to parent or already seen
                if gen_canon == parent_canon or gen_canon in seen:
                    continue
                seen.add(gen_canon)

                gen_feats = count_features(gen_mol)
                # Check that the target descriptor changed in the intended direction.
                # When strict=True, also verify that no other descriptor has changed.
                other_unchanged = all(
                    gen_feats[i] == parent_feats[i]
                    for i in range(len(parent_feats)) if i != feat_idx
                )
                if (direction * gen_feats[feat_idx] > direction * parent_feats[feat_idx]
                        and (not strict or other_unchanged)):
                    results.append({
                        "mode": mode,
                        "site": site,
                        "fragment": frag_smiles,
                        "parent_smiles": parent_canon,
                        "generated_smiles": gen_canon,
                        "parent_features": parent_feats,
                        "generated_features": gen_feats,
                    })
            except Exception:
                continue

            if len(results) >= max_candidates:
                return results

    return results


# ──────────────────────────────────────────────
# 7. Display results
# ──────────────────────────────────────────────
FEATURE_NAMES = ["NumHAcceptors", "NumHDonors", "NumAromaticRings", "NumRotatableBonds", "NumHeteroatoms"]
# Corresponds to the order of values returned by count_features()

def print_proposals(results):
    if not results:
        print("No valid candidates found.")
        return

    for i, rec in enumerate(results, start=1):
        print(f"=== Candidate {i} ===")
        print(f"mode: {rec['mode']},  site: {rec['site']},  fragment: {rec['fragment']}")
        print(f"parent SMILES   : {rec['parent_smiles']}")
        print(f"generated SMILES: {rec['generated_smiles']}")
        print(f"features {FEATURE_NAMES}")
        print(f"  parent   : {rec['parent_features']}")
        print(f"  generated: {rec['generated_features']}")
        print()
