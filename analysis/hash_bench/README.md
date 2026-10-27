# Structural hash benchmark

Standalone analysis of algorithms for **constitutional-isomer identity**.
Nothing in `src/` is imported or modified; the production ReaxTools identity
model is left untouched.

Question: which algorithm can produce a *true* structural identifier, and what
does it cost?

## Quick start

```bash
python3 analysis/hash_bench/run_benchmark.py
```

Outputs are written to `analysis/hash_bench/results/`:

| file | contents |
| --- | --- |
| `accuracy_groups.csv` | per-group unique-ID counts and pair confusion |
| `accuracy_totals.csv` | realistic / adversarial / all pair accuracy |
| `timing.csv` | best-of-5 wall time per method and graph size |
| `time_scaling.png` | log-log timing plot |
| `permutation_invariance.csv` | atom-order sensitivity |
| `summary.md`, `summary.json` | generated report and machine-readable data |

## What is compared

All methods operate on the same **element-labelled heavy-atom graph with a
per-atom hydrogen count**. A hydrogen is a leaf, so this is exactly the
constitutional graph; the representation removes the arbitrary choice of
explicit versus implicit H without losing information.

| method | invariant / canonical form |
| --- | --- |
| `formula` | molecular formula (the current ReaxTools model) |
| `wl1` | 1-dimensional Weisfeiler-Lehman / Morgan colour refinement |
| `wl2` | 2-dimensional WL on ordered vertex pairs |
| `canonical_smiles` | RDKit canonical SMILES, non-isomeric |
| `inchi` | standard InChI |
| `brute_force` | lexicographically minimal labelled adjacency string (audit only) |

Ground truth is **VF2 graph isomorphism** on the labelled heavy-atom graph.
For groups with at most eight heavy atoms, a complete brute-force canonical
form is used as an independent audit; RDKit canonical SMILES, VF2 and the
brute-force form agree on every such group.

## Test sets

1. **Curated constitutional isomers** - C2H6O, C3H6O, C3H8O, C4H10O, C4H8,
   C4H6, C5H12, C6H6, C5H5NO. Includes ethanol/dimethyl ether and the
   2-pyridone/2-hydroxypyridine tautomer pair.
2. **Complete alkane sets** - all acyclic constitutional isomers of C7H16,
   C8H18, C9H20 and C10H22 (9, 18, 35 and 75 structures), generated from all
   non-isomorphic trees of maximum degree four.
3. **Adversarial regular graphs** - every sampled isomorphism class of
   connected 3-regular C_nH_n and 4-regular C_n skeletons for n = 6-12. These
   are valence-valid molecules, but deliberately highly symmetric.
4. **Strongly regular pair** - Shrikhande graph and the 4x4 rook graph, built
   as 6-regular S16 molecules. They are non-isomorphic but have identical
   2-WL colour histograms.

## Accuracy

Pair accuracy treats two molecules as the same only when they are isomorphic
as labelled heavy-atom graphs. Every pair here is same-formula by
construction, so the formula baseline scores 0: this metric measures isomer
discrimination, not general species-counting accuracy.

### Realistic organic sets (curated + alkanes + tautomer pair)

| algorithm | pairs | false merges | false splits | pair accuracy |
| --- | ---: | ---: | ---: | ---: |
| formula | 3644 | 3644 | 0 | 0.0000 |
| 1-WL | 3644 | 1 | 0 | 0.9997 |
| 2-WL | 3644 | 0 | 0 | 1.0000 |
| canonical SMILES | 3644 | 0 | 0 | 1.0000 |
| InChI | 3644 | 1 | 0 | 0.9997 |

The single 1-WL failure is the C6H6 pair K(3,3) versus the triangular prism:
both carbon skeletons are 3-regular on six vertices, so colour refinement
never separates them. 2-WL separates this pair.

The single InChI failure is 2-pyridone versus 2-hydroxypyridine: standard
InChI applies mobile-H/tautomer normalisation and returns the same string.

### Adversarial regular graphs

| algorithm | pairs | false merges | false splits | pair accuracy |
| --- | ---: | ---: | ---: | ---: |
| formula | 4614 | 4614 | 0 | 0.0000 |
| 1-WL | 4614 | 4614 | 0 | 0.0000 |
| 2-WL | 4614 | 1 | 0 | 0.9998 |
| canonical SMILES | 4614 | 0 | 0 | 1.0000 |
| InChI | 4614 | 0 | 0 | 1.0000 |

1-WL cannot distinguish any two non-isomorphic d-regular graphs with the same
n and d. The single 2-WL failure is the Shrikhande/rook pair: same strongly
regular parameters, identical pair-colour histogram. Canonical labelling
separates it.

### All pairs

| algorithm | pairs | false merges | false splits | pair accuracy |
| --- | ---: | ---: | ---: | ---: |
| formula | 8258 | 8258 | 0 | 0.0000 |
| 1-WL | 8258 | 4615 | 0 | 0.4411 |
| 2-WL | 8258 | 1 | 0 | 0.9999 |
| canonical SMILES | 8258 | 0 | 0 | 1.0000 |
| InChI | 8258 | 1 | 0 | 0.9999 |

The all-pairs 1-WL number is dominated by the adversarial regular sets; on
the realistic sets 1-WL is 99.97% accurate. That is the central practical
trade-off: 1-WL is cheap and almost always right for valence-4 organic
molecules, but it is not exact and has a real C6H6 counterexample.

## Timing

Best of five calls, Python 3.11 + RDKit 2022.09 on this machine. Times are
per molecule and exclude molecule construction/parsing.

### Linear alkane chain

| heavy atoms | formula | 1-WL full (rounds) | 1-WL k=4 | canonical SMILES | InChI |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 0.001 ms | 0.060 ms (3) | 0.060 ms | 0.008 ms | 0.061 ms |
| 32 | 0.002 ms | 0.631 ms (15) | 0.222 ms | 0.042 ms | 0.111 ms |
| 128 | 0.004 ms | 9.58 ms (63) | 0.809 ms | 0.332 ms | 0.430 ms |
| 512 | 0.012 ms | 147.7 ms (255) | 3.39 ms | 3.82 ms | 3.84 ms |

Empirical log-log slopes: 1-WL full **1.89**, 1-WL k=4 **0.95**,
canonical SMILES **1.49**, InChI **0.99**.

Full 1-WL is quadratic on chains because colours propagate one bond per round
and a path has diameter n. Capping at four rounds restores linear cost, but
four rounds are a weaker invariant (they still cannot separate regular
graphs). Canonical SMILES is as fast as k=4 WL at 512 atoms and is exact.

### Random 3-regular cage

| heavy atoms | 1-WL full (rounds) | 1-WL k=4 | canonical SMILES | InChI | 2-WL |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 0.044 ms (1) | 0.043 ms | 0.015 ms | 0.090 ms | 0.78 ms |
| 16 | 0.074 ms (1) | 0.070 ms | 0.032 ms | 0.114 ms | 7.62 ms |
| 24 | 0.108 ms (1) | 0.094 ms | 0.058 ms | 0.168 ms | 26.4 ms |
| 32 | 0.117 ms (1) | 0.121 ms | 0.081 ms | 0.223 ms | - |

2-WL empirical slope on these graphs is **3.19** (O(n^3) per round, constant
number of rounds), so it is not practical beyond roughly 50-100 heavy atoms
in Python. 1-WL stabilises in one round here because all vertices start with
the same colour - which is also why it is useless on these graphs.

## Why a fixed-length hash is not automatically "true"

A hash is only as good as its input. The digest length controls accidental
collisions, not structural resolution:

| species count | 32-bit | 64-bit | 128-bit |
| ---: | ---: | ---: | ---: |
| 10^3 | 1.2e-4 | 2.7e-14 | ~0 |
| 10^4 | 1.2e-2 | 2.7e-12 | ~0 |
| 10^5 | 6.9e-1 | 2.7e-10 | ~0 |
| 10^6 | ~1 | 2.7e-8 | ~0 |
| 10^7 | ~1 | 2.7e-6 | ~0 |

(p = 1 - exp(-n(n-1)/2^(b+1))). The current ReaxTools ID is 32-bit and its
equality test compares only that integer, so collision means silent species
merging. Even after switching to a structural invariant, use at least 128
bits (or a 256-bit digest truncated to 128).

## Recommendation

1. **Do not use 1-WL alone as the final identity.** It is exact on almost all
   organic molecules, but K(3,3) versus prismane is a valence-valid C6H6
   counterexample, and every d-regular skeleton collides by construction.
2. **Use a canonical labelling as ground truth.** Practical options:
   - RDKit canonical SMILES (non-isomeric) - easiest, proven, but adds a
     large dependency to the C++ core;
   - nauty/bliss-style individualisation-refinement in C++ - exact,
     dependency-free, empirically near-linear on molecular graphs;
   - hybrid: use 1-WL/2-WL only as a bucket key, then verify with exact graph
     isomorphism (VF2) against the bucket representatives, or canonicalise
     every new candidate. WL is a filter, not an identity: do not return the
     WL hash as the final ID without an exact fallback.
3. **Hash the canonical form with at least 128 bits.** Keep the formula for
   display, store the structural digest as the identity.
4. **Exclude bond order** unless it is canonicalised. ReaxTools assigns bond
   orders by a greedy pass over an unordered bond set; including them would
   make the identity depend on insertion order. Constitutional isomerism is
   already captured by connectivity plus element/H labels.
5. **Do not use standard InChI as an exact constitutional identifier.**
   Mobile-H normalisation merges tautomers (2-pyridone/2-hydroxypyridine).
   Use non-isomeric canonical SMILES, or fixed-H InChI.
6. **Stereoisomers are a separate switch.** Non-isomeric canonical SMILES
   correctly merges cis/trans and R/S for constitutional identity; use the
   isomeric form if stereo matters.

## Caveats

- The conventional organic sets are small; the regular-graph sets are
  deliberately adversarial, and most are highly strained or hypothetical
  (prismane is real, K(3,3)-C6H6 is not a known isolable molecule). They are
  counterexamples to completeness, not a claim about typical ReaxFF output.
- VF2 ground truth ignores bond order, charge and spin. That is intentional
  for constitutional identity; a production identity model should decide
  explicitly whether charge/spin belong in the key.
- RDKit canonical SMILES is the de facto standard, not a formally proven
  canonical form for all graphs. The brute-force audit covers n <= 8 and the
  VF2 checks cover every test group; no discrepancy was observed.
- Timings are Python/RDKit. A C++ 1-WL implementation would be much faster in
  absolute terms, but the asymptotic shapes (quadratic full WL on chains,
  cubic 2-WL, near-linear canonicalisation on molecules) are the relevant
  result.
- This environment had no external network access, so the complexity
  statements are textbook ones rather than freshly cited references.
