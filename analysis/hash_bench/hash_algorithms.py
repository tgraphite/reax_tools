"""Standalone structural-identity algorithms for molecular graphs.

Everything here operates on an RDKit molecule with *implicit* hydrogens.
The graph used for identity is the heavy-atom graph plus a per-atom H count:
a hydrogen is a leaf, so (element, total-H) on each heavy atom plus the
heavy-atom connectivity is exactly the constitutional graph.

No production code in src/ is imported or modified.
"""
from __future__ import annotations

import hashlib
from collections import Counter
from itertools import permutations

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def graph_of(mol: Chem.Mol):
    """Return (labels, adjacency) for the heavy-atom graph with H counts."""
    n = mol.GetNumAtoms()
    labels = [(atom.GetSymbol(), atom.GetTotalNumHs()) for atom in mol.GetAtoms()]
    adj = [set() for _ in range(n)]
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adj[i].add(j)
        adj[j].add(i)
    return labels, adj


def mol_from_graph(n: int, edges, symbol: str = "C") -> Chem.Mol:
    """Build a sanitized molecule from an unlabeled simple graph."""
    rw = Chem.RWMol()
    for _ in range(n):
        rw.AddAtom(Chem.Atom(symbol))
    for i, j in edges:
        rw.AddBond(int(i), int(j), Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return mol


# --------------------------------------------------------------------------
# baseline: the current ReaxTools identity model
# --------------------------------------------------------------------------
def formula_id(mol: Chem.Mol) -> str:
    """Formula-level identity: all constitutional isomers collide."""
    return sha256_text("formula|" + rdMolDescriptors.CalcMolFormula(mol))


def _stable_color(value) -> int:
    """Graph-independent 64-bit colour id (never a process-local Python hash)."""
    return int.from_bytes(
        hashlib.blake2b(repr(value).encode("utf-8"), digest_size=8).digest(),
        "little",
    )


def _partition_stable(old_colors, new_colors) -> bool:
    """True when no old colour class was split by the update."""
    old_classes = set(old_colors.values()) if isinstance(old_colors, dict) else set(old_colors)
    mapping = (
        set(zip(old_colors.values(), new_colors.values()))
        if isinstance(old_colors, dict)
        else set(zip(old_colors, new_colors))
    )
    return len(mapping) == len(old_classes)


# --------------------------------------------------------------------------
# 1-WL / colour refinement (Morgan-style)
# --------------------------------------------------------------------------
def wl1_hash(mol: Chem.Mol, max_rounds: int | None = None):
    """1-dimensional Weisfeiler-Lehman hash.

    Returns (hex_digest, rounds).  The invariant is the multiset of final
    vertex colours.  Isomorphic graphs always agree; non-isomorphic graphs
    may also agree (e.g. any two d-regular graphs with the same n and d).
    """
    labels, adj = graph_of(mol)
    n = len(labels)
    m = sum(len(a) for a in adj) // 2

    color = [_stable_color(("init", lab)) for lab in labels]

    rounds = 0
    limit = n if max_rounds is None else max_rounds
    for rounds in range(1, limit + 1):
        signature = [
            (color[i], tuple(sorted(color[j] for j in adj[i]))) for i in range(n)
        ]
        new_color = [_stable_color(sig) for sig in signature]
        stable = _partition_stable(color, new_color)
        color = new_color
        if stable:
            break

    histogram = tuple(sorted(Counter(color).items()))
    digest = sha256_text(f"wl1|n={n}|m={m}|{histogram}")
    return digest, rounds


# --------------------------------------------------------------------------
# 2-WL (pair refinement)
# --------------------------------------------------------------------------
def wl2_hash(mol: Chem.Mol, max_rounds: int | None = None):
    """2-dimensional Weisfeiler-Lehman hash on ordered vertex pairs.

    Cost is O(rounds * n^3 log n) time and O(n^2) memory, so this is only
    benchmarked for small graphs.  2-WL is strictly stronger than 1-WL but
    is still not a complete invariant (e.g. pairs of strongly regular graphs
    with the same parameters).
    """
    labels, adj = graph_of(mol)
    n = len(labels)
    adjacent = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in adj[i]:
            adjacent[i][j] = 1
    m = sum(sum(row) for row in adjacent) // 2

    color = {
        (u, v): _stable_color(("init", labels[u], labels[v], adjacent[u][v]))
        for u in range(n)
        for v in range(n)
    }

    rounds = 0
    limit = n if max_rounds is None else max_rounds
    for rounds in range(1, limit + 1):
        signature = {}
        for u in range(n):
            for v in range(n):
                multiset = tuple(
                    sorted((color[(u, w)], color[(w, v)]) for w in range(n))
                )
                signature[(u, v)] = (color[(u, v)], multiset)
        new_color = {pair: _stable_color(sig) for pair, sig in signature.items()}
        stable = _partition_stable(color, new_color)
        color = new_color
        if stable:
            break

    histogram = tuple(sorted(Counter(color.values()).items()))
    digest = sha256_text(f"wl2|n={n}|m={m}|{histogram}")
    return digest, rounds


# --------------------------------------------------------------------------
# canonical representations
# --------------------------------------------------------------------------
def canonical_smiles_id(mol: Chem.Mol) -> str:
    """RDKit canonical SMILES without stereochemistry.

    This is the practical reference implementation of a canonical form for
    constitutional identity.  It is not a fixed-length digest; callers can
    hash the returned string if a compact ID is needed.
    """
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)


def inchi_id(mol: Chem.Mol) -> str:
    """Standard InChI.

    Note: standard InChI applies mobile-H/tautomer normalisation, so some
    constitutional isomers (e.g. 2-pyridone and 2-hydroxypyridine) receive
    the same InChI.  It is therefore not exact for constitutional identity.
    """
    value = Chem.MolToInchi(mol)
    return value if value else "<inchi-failed>"


# --------------------------------------------------------------------------
# brute-force canonical form (ground truth for tiny graphs)
# --------------------------------------------------------------------------
def brute_force_id(mol: Chem.Mol) -> str:
    """Lexicographically minimal labelled adjacency string.

    Complete canonical form, O(n! * n^2).  Only usable for n <= 9 heavy
    atoms; it is used to audit the other canonicalisers.
    """
    labels, adj = graph_of(mol)
    n = len(labels)
    if n > 9:
        raise ValueError("brute_force_id is limited to 9 heavy atoms")

    best = None
    for perm in permutations(range(n)):
        text = tuple(labels[perm[i]] for i in range(n)) + tuple(
            1 if perm[j] in adj[perm[i]] else 0
            for i in range(n)
            for j in range(i + 1, n)
        )
        if best is None or text < best:
            best = text
    return repr(best)


# --------------------------------------------------------------------------
# registry used by the benchmark driver
# --------------------------------------------------------------------------
def id_functions(include_wl2: bool = True):
    funcs = {
        "formula": lambda mol: formula_id(mol),
        "wl1": lambda mol: wl1_hash(mol)[0],
        "canonical_smiles": canonical_smiles_id,
        "inchi": inchi_id,
    }
    if include_wl2:
        funcs["wl2"] = lambda mol: wl2_hash(mol)[0]
    return funcs