"""Benchmark structural-identity algorithms on constitutional isomers.

Run:
    python3 analysis/hash_bench/run_benchmark.py

Outputs (analysis/hash_bench/results/):
    accuracy_groups.csv   per-group unique-ID counts and pair confusion
    accuracy_totals.csv   aggregate pair confusion per algorithm
    timing.csv            median wall time per algorithm and graph size
    time_scaling.png      log-log timing plot
    summary.json          machine-readable summary
    summary.md            human-readable summary

Nothing in src/ is imported or modified.
"""
from __future__ import annotations

import csv
import json
import statistics
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from datasets import (  # noqa: E402
    alkane_families,
    curated_families,
    regular_classes,
    strong_regular_pair,
    tautomer_pair,
)
from hash_algorithms import (  # noqa: E402
    brute_force_id,
    canonical_smiles_id,
    formula_id,
    inchi_id,
    mol_from_graph,
    wl1_hash,
    wl2_hash,
)

RESULTS = HERE / "results"
RESULTS.mkdir(exist_ok=True)

ALGOS = ["formula", "wl1", "wl2", "canonical_smiles", "inchi"]


# --------------------------------------------------------------------------
# graph helpers
# --------------------------------------------------------------------------
def to_nx(mol: Chem.Mol) -> nx.Graph:
    graph = nx.Graph()
    for atom in mol.GetAtoms():
        graph.add_node(
            atom.GetIdx(), element=atom.GetSymbol(), h=atom.GetTotalNumHs()
        )
    for bond in mol.GetBonds():
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    return graph


def node_match(a, b):
    return a["element"] == b["element"] and a["h"] == b["h"]


def truth_partition(mols):
    """VF2-based isomorphism classes; independent of all tested algorithms."""
    n = len(mols)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    graphs = [to_nx(mol) for mol in mols]
    for i in range(n):
        for j in range(i + 1, n):
            if find(i) == find(j):
                continue
            if nx.is_isomorphic(graphs[i], graphs[j], node_match=node_match):
                union(i, j)
    return [find(i) for i in range(n)]


def ids_for(mols, include_wl2=True):
    out = {
        "formula": [formula_id(mol) for mol in mols],
        "wl1": [wl1_hash(mol)[0] for mol in mols],
        "canonical_smiles": [canonical_smiles_id(mol) for mol in mols],
        "inchi": [inchi_id(mol) for mol in mols],
    }
    if include_wl2:
        out["wl2"] = [wl2_hash(mol)[0] for mol in mols]
    return out


# --------------------------------------------------------------------------
# accuracy
# --------------------------------------------------------------------------
def evaluate_group(name, mols, include_wl2=True, category="realistic"):
    mols = [mol for mol in mols if mol is not None]
    if len(mols) < 2:
        return None

    truth = truth_partition(mols)
    n_struct = len(set(truth))
    ids = ids_for(mols, include_wl2=include_wl2)

    row = {
        "group": name,
        "category": category,
        "n_molecules": len(mols),
        "n_structures": n_struct,
        "pairs": len(mols) * (len(mols) - 1) // 2,
    }
    for algo in ALGOS:
        if algo not in ids:
            continue
        values = ids[algo]
        row[f"n_ids_{algo}"] = len(set(values))

    # pair confusion
    for algo in ALGOS:
        if algo not in ids:
            continue
        false_merge = false_split = 0
        for i in range(len(mols)):
            for j in range(i + 1, len(mols)):
                truly_same = truth[i] == truth[j]
                predicted_same = ids[algo][i] == ids[algo][j]
                if truly_same and not predicted_same:
                    false_split += 1
                elif (not truly_same) and predicted_same:
                    false_merge += 1
        row[f"fm_{algo}"] = false_merge
        row[f"fs_{algo}"] = false_split
    return row


def accuracy_section():
    groups = {}
    groups.update(curated_families())
    groups.update(alkane_families([7, 8, 9, 10]))

    # Adversarial 1-WL sets: connected d-regular carbon skeletons.
    for n, d in [(6, 3), (8, 3), (10, 3), (12, 3), (8, 4), (10, 4)]:
        groups[f"regular_d{d}_n{n}"] = regular_classes(
            n, d, seeds=range(0, 500)
        )

    # 2-WL counterexample (valid S16 molecules).
    try:
        sr = strong_regular_pair()
        if all(mol is not None for mol in sr):
            groups["strong_regular_S16"] = sr
    except Exception as exc:  # pragma: no cover - environment dependent
        print(f"[warn] strong-regular pair skipped: {exc}")

    rows = []
    for name, mols in groups.items():
        category = (
            "adversarial"
            if name.startswith("regular_") or name.startswith("strong_regular")
            else "realistic"
        )
        row = evaluate_group(name, mols, category=category)
        if row:
            rows.append(row)

    with open(RESULTS / "accuracy_groups.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    totals = {}
    for category in ["realistic", "adversarial", "all"]:
        totals[category] = {
            algo: {"false_merge": 0, "false_split": 0, "pairs": 0}
            for algo in ALGOS
        }
    for row in rows:
        categories = [row["category"], "all"]
        for category in categories:
            for algo in ALGOS:
                if f"fm_{algo}" in row:
                    totals[category][algo]["false_merge"] += row[f"fm_{algo}"]
                    totals[category][algo]["false_split"] += row[f"fs_{algo}"]
                    totals[category][algo]["pairs"] += row["pairs"]

    with open(RESULTS / "accuracy_totals.csv", "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "category",
                "algorithm",
                "pairs",
                "false_merge",
                "false_split",
                "pair_accuracy",
            ],
        )
        writer.writeheader()
        for category in ["realistic", "adversarial", "all"]:
            for algo in ALGOS:
                t = totals[category][algo]
                if t["pairs"] == 0:
                    continue
                accuracy = 1.0 - (t["false_merge"] + t["false_split"]) / t["pairs"]
                writer.writerow(
                    {
                        "category": category,
                        "algorithm": algo,
                        "pairs": t["pairs"],
                        "false_merge": t["false_merge"],
                        "false_split": t["false_split"],
                        "pair_accuracy": f"{accuracy:.6f}",
                    }
                )

    return rows, totals


def brute_force_audit():
    """Check canonical SMILES / VF2 against a complete canonical form."""
    groups = {}
    groups.update(curated_families())
    groups["regular_d3_n6"] = regular_classes(6, 3, range(0, 50))
    groups["regular_d3_n8"] = regular_classes(8, 3, range(0, 200))

    report = []
    for name, mols in groups.items():
        mols = [mol for mol in mols if mol is not None]
        if len(mols) < 2 or any(mol.GetNumAtoms() > 8 for mol in mols):
            continue
        bf = [brute_force_id(mol) for mol in mols]
        cs = [canonical_smiles_id(mol) for mol in mols]
        vf2 = truth_partition(mols)
        report.append(
            {
                "group": name,
                "n_molecules": len(mols),
                "n_bruteforce": len(set(bf)),
                "n_smiles": len(set(cs)),
                "n_vf2": len(set(vf2)),
                "agree": len(set(bf)) == len(set(cs)) == len(set(vf2)),
            }
        )
    return report


def permutation_invariance(trials=5):
    """Check that each identity is invariant to input atom ordering."""
    import random

    random.seed(20260804)
    groups = {}
    groups.update(curated_families())
    groups.update(alkane_families([8]))
    groups["regular_d3_n8"] = regular_classes(8, 3, range(0, 100))

    sample = []
    for mols in groups.values():
        sample.extend(mols[:4])

    funcs = {
        "formula": lambda mol: formula_id(mol),
        "wl1": lambda mol: wl1_hash(mol)[0],
        "wl2": lambda mol: wl2_hash(mol)[0],
        "canonical_smiles": canonical_smiles_id,
        "inchi": inchi_id,
    }
    report = {algo: {"molecules": 0, "variants": 0, "false_split": 0} for algo in funcs}

    for mol in sample:
        if mol is None:
            continue
        n = mol.GetNumAtoms()
        variants = [mol]
        for _ in range(trials):
            order = list(range(n))
            random.shuffle(order)
            variants.append(Chem.RenumberAtoms(mol, order))

        for algo, fn in funcs.items():
            values = {fn(variant) for variant in variants}
            report[algo]["molecules"] += 1
            report[algo]["variants"] += len(variants)
            if len(values) > 1:
                report[algo]["false_split"] += 1

    with open(RESULTS / "permutation_invariance.csv", "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["algorithm", "molecules", "variants", "false_split"],
        )
        writer.writeheader()
        for algo, item in report.items():
            writer.writerow({"algorithm": algo, **item})
    return report


# --------------------------------------------------------------------------
# timing
# --------------------------------------------------------------------------
def chain_mol(n_heavy: int) -> Chem.Mol:
    mol = Chem.MolFromSmiles("C" * n_heavy)
    if mol is None:
        raise RuntimeError(f"failed to build C{n_heavy} chain")
    return mol


def regular_mol(n: int, d: int = 3) -> Chem.Mol:
    for seed in range(1000):
        graph = nx.random_regular_graph(d, n, seed=seed)
        if nx.is_connected(graph):
            return mol_from_graph(n, graph.edges())
    raise RuntimeError("could not build a connected regular graph")


def timed(fn, mol, repeat=5, inner=None):
    """Return the best per-call wall time (seconds).

    Small molecules are timed in batches so that sub-microsecond calls are
    not dominated by clock resolution; the minimum over repeats is used to
    suppress scheduler noise.
    """
    fn(mol)  # warm-up
    if inner is None:
        inner = max(1, min(200, 2000 // max(mol.GetNumAtoms(), 1)))
    samples = []
    for _ in range(repeat):
        start = time.perf_counter()
        for _ in range(inner):
            fn(mol)
        samples.append((time.perf_counter() - start) / inner)
    return min(samples)


def timing_section():
    rows = []

    configurations = []
    for n in [8, 16, 32, 64, 128, 256, 512]:
        configurations.append(("chain", n, chain_mol(n)))
    for n in [8, 12, 16, 24, 32]:
        configurations.append(("regular3", n, regular_mol(n, 3)))

    for family, n_heavy, mol in configurations:
        n_total = mol.GetNumAtoms()
        base = {
            "family": family,
            "n_heavy": n_heavy,
            "n_total": n_total,
        }

        rows.append({**base, "method": "formula", "best_ms": timed(formula_id, mol) * 1e3, "rounds": ""})

        digest, rounds = wl1_hash(mol)
        rows.append({**base, "method": "wl1_full", "best_ms": timed(lambda m: wl1_hash(m)[0], mol) * 1e3, "rounds": rounds})
        rows.append({**base, "method": "wl1_k4", "best_ms": timed(lambda m: wl1_hash(m, max_rounds=4)[0], mol) * 1e3, "rounds": 4})

        rows.append({**base, "method": "canonical_smiles", "best_ms": timed(canonical_smiles_id, mol) * 1e3, "rounds": ""})
        rows.append({**base, "method": "inchi", "best_ms": timed(inchi_id, mol) * 1e3, "rounds": ""})

        if n_heavy <= 24:
            _, rounds2 = wl2_hash(mol)
            rows.append(
                {
                    **base,
                    "method": "wl2_full",
                    "best_ms": timed(lambda m: wl2_hash(m)[0], mol, inner=1) * 1e3,
                    "rounds": rounds2,
                }
            )

    with open(RESULTS / "timing.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def plot_timing(rows):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for ax, family, title in zip(
        axes, ["chain", "regular3"], ["linear alkane chain", "random 3-regular cage"]
    ):
        for method in ["wl1_full", "wl1_k4", "canonical_smiles", "inchi", "wl2_full", "formula"]:
            subset = [r for r in rows if r["family"] == family and r["method"] == method]
            if not subset:
                continue
            xs = [r["n_heavy"] for r in subset]
            ys = [max(r["best_ms"], 1e-4) for r in subset]
            ax.plot(xs, ys, marker="o", linewidth=1.4, label=method)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("heavy atoms")
        ax.set_title(title)
        ax.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel("best wall time (ms, best of 5)")
    axes[0].legend(fontsize=8)
    fig.suptitle("Structural identity cost per molecule (log-log)")
    fig.tight_layout()
    fig.savefig(RESULTS / "time_scaling.png", dpi=180)
    plt.close(fig)


# --------------------------------------------------------------------------
# summary
# --------------------------------------------------------------------------
def scaling_exponents(rows):
    out = []
    for family in ["chain", "regular3"]:
        for method in ["wl1_full", "wl1_k4", "canonical_smiles", "inchi", "wl2_full"]:
            subset = [
                r
                for r in rows
                if r["family"] == family
                and r["method"] == method
                and r["best_ms"] > 0
            ]
            if len(subset) < 3:
                continue
            xs = np.log([r["n_heavy"] for r in subset])
            ys = np.log([r["best_ms"] for r in subset])
            slope = float(np.polyfit(xs, ys, 1)[0])
            out.append({"family": family, "method": method, "loglog_slope": round(slope, 3)})
    return out


def write_summary(acc_rows, totals, timing_rows, audit, perm_report, tautomer_note):
    exponents = scaling_exponents(timing_rows)

    lines = []
    lines.append("# Structural hash benchmark\n")
    lines.append(
        "Standalone benchmark for constitutional-isomer identity. "
        "No production code is imported or modified.\n"
    )

    lines.append("## Accuracy on constitutional-isomer groups\n")
    lines.append(
        "Pair accuracy treats two molecules as the same only when they are "
        "isomorphic as labelled heavy-atom graphs (VF2 ground truth).\n"
    )
    for category in ["realistic", "adversarial", "all"]:
        lines.append(f"\n### {category} pairs\n")
        lines.append("| algorithm | pairs | false merges | false splits | pair accuracy |")
        lines.append("| --- | ---: | ---: | ---: | ---: |")
        for algo in ALGOS:
            t = totals[category][algo]
            if t["pairs"] == 0:
                continue
            accuracy = 1.0 - (t["false_merge"] + t["false_split"]) / t["pairs"]
            lines.append(
                f"| {algo} | {t['pairs']} | {t['false_merge']} | {t['false_split']} | {accuracy:.4f} |"
            )

    lines.append("\n## Per-group detail\n")
    lines.append("| group | structures | formula | wl1 | wl2 | canonical SMILES | InChI |")
    lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for row in acc_rows:
        lines.append(
            "| {group} | {n_structures} | {formula} | {wl1} | {wl2} | {cs} | {inchi} |".format(
                group=row["group"],
                n_structures=row["n_structures"],
                formula=row.get("n_ids_formula", ""),
                wl1=row.get("n_ids_wl1", ""),
                wl2=row.get("n_ids_wl2", ""),
                cs=row.get("n_ids_canonical_smiles", ""),
                inchi=row.get("n_ids_inchi", ""),
            )
        )

    lines.append("\n## Atom-order invariance\n")
    lines.append("Each molecule was renumbered five times; a non-zero count means the ID depends on input atom order.\n")
    lines.append("| algorithm | molecules | variants | molecules with false splits |")
    lines.append("| --- | ---: | ---: | ---: |")
    for algo, item in perm_report.items():
        lines.append(
            f"| {algo} | {item['molecules']} | {item['variants']} | {item['false_split']} |"
        )

    lines.append("\n## Timing scaling\n")
    lines.append("| family | method | log-log slope |")
    lines.append("| --- | --- | ---: |")
    for item in exponents:
        lines.append(f"| {item['family']} | {item['method']} | {item['loglog_slope']} |")

    lines.append("\n## Brute-force audit (n <= 8 heavy atoms)\n")
    lines.append("| group | molecules | brute force classes | canonical SMILES classes | VF2 classes | agree |")
    lines.append("| --- | ---: | ---: | ---: | ---: | --- |")
    for item in audit:
        lines.append(
            f"| {item['group']} | {item['n_molecules']} | {item['n_bruteforce']} | "
            f"{item['n_smiles']} | {item['n_vf2']} | {item['agree']} |"
        )

    lines.append("\n## Notes\n")
    lines.append(tautomer_note)
    lines.append(
        "\n1-WL is exact on trees and on all non-regular graphs in these sets, "
        "but it is provably incomplete: every pair of non-isomorphic d-regular "
        "graphs with the same n and d has the same 1-WL hash. 2-WL fixes the "
        "6-vertex K(3,3)/prism case but still collides on the Shrikhande and "
        "4x4 rook graphs (6-regular S16). Only canonical labelling is exact.\n"
    )

    (RESULTS / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    (RESULTS / "summary.json").write_text(
        json.dumps(
            {
                "totals": totals,
                "scaling": exponents,
                "brute_force_audit": audit,
                "permutation_invariance": perm_report,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return "\n".join(lines)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main():
    print("=" * 72)
    print("Accuracy on constitutional-isomer groups")
    print("=" * 72)
    acc_rows, totals = accuracy_section()
    for category in ["realistic", "adversarial", "all"]:
        print(f"\n-- {category} pairs --")
        for algo in ALGOS:
            t = totals[category][algo]
            if not t["pairs"]:
                continue
            accuracy = 1.0 - (t["false_merge"] + t["false_split"]) / t["pairs"]
            print(
                f"  {algo:18s} pairs={t['pairs']:6d} "
                f"false_merge={t['false_merge']:6d} false_split={t['false_split']:6d} "
                f"acc={accuracy:.4f}"
            )

    print("\n" + "=" * 72)
    print("Brute-force audit")
    print("=" * 72)
    audit = brute_force_audit()
    for item in audit:
        print(
            f"  {item['group']:24s} bf={item['n_bruteforce']:3d} "
            f"smiles={item['n_smiles']:3d} vf2={item['n_vf2']:3d} "
            f"agree={item['agree']}"
        )

    print("\n" + "=" * 72)
    print("Atom-order invariance")
    print("=" * 72)
    perm_report = permutation_invariance()
    for algo, item in perm_report.items():
        print(
            f"  {algo:18s} molecules={item['molecules']:3d} "
            f"variants={item['variants']:3d} false_split={item['false_split']}"
        )

    print("\n" + "=" * 72)
    print("Timing")
    print("=" * 72)
    timing_rows = timing_section()
    plot_timing(timing_rows)

    tautomers = tautomer_pair()
    note = (
        "Standard InChI merges 2-pyridone and 2-hydroxypyridine "
        f"({'same' if inchi_id(tautomers[0]) == inchi_id(tautomers[1]) else 'different'} InChI). "
        "Use non-isomeric canonical SMILES, or fixed-H InChI, when "
        "constitutional identity must be exact."
    )
    print("\n" + note)

    write_summary(acc_rows, totals, timing_rows, audit, perm_report, note)
    print(f"\nResults written to {RESULTS}")


if __name__ == "__main__":
    main()
