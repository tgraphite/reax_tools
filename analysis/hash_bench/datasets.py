"""Test-set construction for the structural-hash benchmark.

All sets are groups of constitutional isomers (same molecular formula,
different heavy-atom connectivity / H placement).  Ground truth is graph
isomorphism on the labelled heavy-atom graph (element, total H), never the
output of the algorithm being tested.
"""
from __future__ import annotations

import networkx as nx
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from hash_algorithms import mol_from_graph


CURATED = {
    "C2H6O": ["CCO", "COC"],
    "C3H6O": [
        "CCC=O",          # propanal
        "CC(C)=O",        # acetone
        "C=CCO",          # allyl alcohol
        "CC1CO1",         # propylene oxide
        "COC=C",          # methyl vinyl ether
        "C1CC1O",         # cyclopropanol
    ],
    "C3H8O": ["CCCO", "CC(O)C", "CCOC"],
    "C4H10O": [
        "CCCCO",
        "CC(O)CC",
        "CC(C)CO",
        "CC(C)(C)O",
        "CCCOC",
        "CC(C)OC",
        "CCOCC",
    ],
    "C4H8": ["C=CCC", "CC=CC", "CC(C)=C", "C1CCC1", "CC1CC1"],
    "C4H6": [
        "C=CC=C",
        "C#CCC",
        "CC#CC",
        "C1=CCC1",
        "CC1=CC1",
        "C=C1CC1",
        "C1C2C1C2",
    ],
    "C5H12": ["CCCCC", "CC(C)CC", "CC(C)(C)C"],
    "C6H6": [
        "c1ccccc1",       # benzene
        "C=C1C=CC=C1",    # fulvene
        "C1=CC2C=CC12",   # Dewar benzene
    ],
    "C5H5NO": [
        "O=c1cccc[nH]1",  # 2-pyridone
        "Oc1ccccn1",      # 2-hydroxypyridine (standard InChI merges these)
    ],
}


def canonical_key(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)


def dedup(mols):
    seen = {}
    for mol in mols:
        if mol is None:
            continue
        seen.setdefault(canonical_key(mol), mol)
    return list(seen.values())


def curated_families():
    """Formula -> list of distinct constitutional isomers."""
    families = {}
    for expected, smiles_list in CURATED.items():
        mols = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            formula = rdMolDescriptors.CalcMolFormula(mol)
            if formula != expected:
                continue
            mols.append(mol)
        if mols:
            families[expected] = dedup(mols)

    # Add the two 3-regular C6H6 carbon skeletons.  K(3,3) and the
    # triangular prism are non-isomorphic but 1-WL-indistinguishable.
    k33 = mol_from_graph(6, nx.complete_bipartite_graph(3, 3).edges())
    prism = mol_from_graph(6, nx.circular_ladder_graph(3).edges())
    families["C6H6"] = dedup(families.get("C6H6", []) + [k33, prism])
    return families


def alkane_families(sizes):
    """Complete sets of acyclic alkane constitutional isomers.

    NetworkX enumerates all non-isomorphic trees; trees of maximum degree 4
    are exactly the carbon skeletons of C_n H_(2n+2) alkanes.
    """
    families = {}
    for n in sizes:
        mols = []
        for tree in nx.nonisomorphic_trees(n):
            if max(dict(tree.degree()).values()) > 4:
                continue
            mols.append(mol_from_graph(n, tree.edges()))
        if mols:
            formula = rdMolDescriptors.CalcMolFormula(mols[0])
            families[formula] = dedup(mols)
    return families


def regular_classes(n: int, d: int, seeds) -> list[Chem.Mol]:
    """One representative per isomorphism class of connected d-regular graphs.

    Every vertex is carbon, so the formula is C_n H_(n*(4-d)).  These sets
    are the adversarial case for 1-WL: all graphs with the same (n, d) have
    exactly the same 1-WL colour histogram.
    """
    reps = []
    for seed in seeds:
        try:
            graph = nx.random_regular_graph(d, n, seed=seed)
        except nx.NetworkXError:
            continue
        if not nx.is_connected(graph):
            continue
        if any(nx.is_isomorphic(graph, other) for other in reps):
            continue
        reps.append(graph)
    return [mol_from_graph(n, graph.edges()) for graph in reps]


def tautomer_pair():
    """A known standard-InChI tautomer merge, used for the accuracy caveat."""
    return [
        Chem.MolFromSmiles("O=c1cccc[nH]1"),
        Chem.MolFromSmiles("Oc1ccccn1"),
    ]


def strong_regular_pair():
    """Shrikhande vs 4x4 rook graph: 6-regular, 2-WL-indistinguishable.

    Sulfur has valence 6, so these are valid S16 molecules, not just graphs.
    """
    rook = nx.cartesian_product(nx.complete_graph(4), nx.complete_graph(4))
    shrikhande = nx.Graph()
    for u in range(4):
        for v in range(4):
            shrikhande.add_node((u, v))
    generators = [(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)]
    for u in shrikhande.nodes():
        for du, dv in generators:
            v = ((u[0] + du) % 4, (u[1] + dv) % 4)
            if v != u:
                shrikhande.add_edge(u, v)
    rook = nx.convert_node_labels_to_integers(rook)
    shrikhande = nx.convert_node_labels_to_integers(shrikhande)
    return [
        mol_from_graph(rook.number_of_nodes(), rook.edges(), symbol="S"),
        mol_from_graph(shrikhande.number_of_nodes(), shrikhande.edges(), symbol="S"),
    ]
