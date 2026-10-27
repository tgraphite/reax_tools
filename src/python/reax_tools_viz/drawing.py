"""Drawing-only bond-order assignment for canonical SMILES.

The structural identity is deliberately connectivity-only: every bond in the
canonical SMILES is single. RDKit therefore interprets unsaturated fragments as
radicals. For display only, this module applies the same residual-valence greedy
idea as the C++ bond builder: raise bond orders while both endpoints still have
valence left. The resulting bond orders are not used for identity, counting, or
event detection.
"""

from __future__ import annotations


# Defaults mirror src/cpp/argparser.cpp ELEMENT_MAX_VALENCIES. Unknown elements
# fall back to 4, matching System::set_types().
MAX_VALENCE = {
    "H": 1, "He": 0, "Li": 1, "Be": 2, "B": 3, "C": 4, "N": 4, "O": 2,
    "F": 1, "Ne": 0, "Na": 1, "Mg": 2, "Al": 3, "Si": 4, "P": 6, "S": 6,
    "Cl": 7, "Ar": 0, "K": 1, "Ca": 2, "Br": 7, "I": 7,
}


def _bond_type(order: int, Chem):
    if order >= 3:
        return Chem.BondType.TRIPLE
    if order == 2:
        return Chem.BondType.DOUBLE
    return Chem.BondType.SINGLE


def _greedy_bond_order_mol(base, Chem):
    """Return a sanitized molecule with drawing-only multiple bonds."""
    mol = Chem.RWMol(base)
    h_counts = {atom.GetIdx(): atom.GetTotalNumHs() for atom in mol.GetAtoms()}

    bonds: list[tuple[int, int]] = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if mol.GetAtomWithIdx(i).GetSymbol() == "H" or mol.GetAtomWithIdx(j).GetSymbol() == "H":
            continue
        bonds.append((i, j))
    bonds.sort()

    residual = {
        atom.GetIdx(): MAX_VALENCE.get(atom.GetSymbol(), 4)
        - (atom.GetDegree() + h_counts[atom.GetIdx()])
        for atom in mol.GetAtoms()
    }
    orders = {bond: 1 for bond in bonds}

    # Same fixed-order, repeated-pass greedy as System::build_molecules().
    while True:
        gain = 0
        for i, j in bonds:
            if residual[i] > 0 and residual[j] > 0:
                orders[(i, j)] += 1
                residual[i] -= 1
                residual[j] -= 1
                gain += 2
        if gain == 0:
            break

    for (i, j), order in orders.items():
        mol.GetBondBetweenAtoms(i, j).SetBondType(_bond_type(order, Chem))

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H":
            continue
        idx = atom.GetIdx()
        valence = h_counts[idx] + sum(
            bond.GetBondTypeAsDouble() for bond in atom.GetBonds()
        )
        charge = 0
        if atom.GetSymbol() == "N" and valence > 3:
            charge = int(round(valence - 3))
        elif atom.GetSymbol() == "O" and valence > 2:
            charge = int(round(valence - 2))
        atom.SetFormalCharge(charge)
        atom.SetNumRadicalElectrons(0)
        atom.SetNoImplicit(True)

    Chem.SanitizeMol(mol)
    return mol.GetMol()


def prepare_drawing_mol(smiles: str):
    """Parse a connectivity-only SMILES and add drawing-only bond orders."""
    from rdkit import Chem

    base = Chem.MolFromSmiles(smiles)
    if base is None:
        return None
    try:
        return _greedy_bond_order_mol(base, Chem)
    except Exception:
        # A malformed or unusual valence should not stop the montage.
        return base

