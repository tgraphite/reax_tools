# Structural hash benchmark

Standalone benchmark for constitutional-isomer identity. No production code is imported or modified.

## Accuracy on constitutional-isomer groups

Pair accuracy treats two molecules as the same only when they are isomorphic as labelled heavy-atom graphs (VF2 ground truth).


### realistic pairs

| algorithm | pairs | false merges | false splits | pair accuracy |
| --- | ---: | ---: | ---: | ---: |
| formula | 3644 | 3644 | 0 | 0.0000 |
| wl1 | 3644 | 1 | 0 | 0.9997 |
| wl2 | 3644 | 0 | 0 | 1.0000 |
| canonical_smiles | 3644 | 0 | 0 | 1.0000 |
| inchi | 3644 | 1 | 0 | 0.9997 |

### adversarial pairs

| algorithm | pairs | false merges | false splits | pair accuracy |
| --- | ---: | ---: | ---: | ---: |
| formula | 4614 | 4614 | 0 | 0.0000 |
| wl1 | 4614 | 4614 | 0 | 0.0000 |
| wl2 | 4614 | 1 | 0 | 0.9998 |
| canonical_smiles | 4614 | 0 | 0 | 1.0000 |
| inchi | 4614 | 0 | 0 | 1.0000 |

### all pairs

| algorithm | pairs | false merges | false splits | pair accuracy |
| --- | ---: | ---: | ---: | ---: |
| formula | 8258 | 8258 | 0 | 0.0000 |
| wl1 | 8258 | 4615 | 0 | 0.4411 |
| wl2 | 8258 | 1 | 0 | 0.9999 |
| canonical_smiles | 8258 | 0 | 0 | 1.0000 |
| inchi | 8258 | 1 | 0 | 0.9999 |

## Per-group detail

| group | structures | formula | wl1 | wl2 | canonical SMILES | InChI |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| C2H6O | 2 | 1 | 2 | 2 | 2 | 2 |
| C3H6O | 6 | 1 | 6 | 6 | 6 | 6 |
| C3H8O | 3 | 1 | 3 | 3 | 3 | 3 |
| C4H10O | 7 | 1 | 7 | 7 | 7 | 7 |
| C4H8 | 5 | 1 | 5 | 5 | 5 | 5 |
| C4H6 | 7 | 1 | 7 | 7 | 7 | 7 |
| C5H12 | 3 | 1 | 3 | 3 | 3 | 3 |
| C6H6 | 5 | 1 | 4 | 5 | 5 | 5 |
| C5H5NO | 2 | 1 | 2 | 2 | 2 | 1 |
| C7H16 | 9 | 1 | 9 | 9 | 9 | 9 |
| C8H18 | 18 | 1 | 18 | 18 | 18 | 18 |
| C9H20 | 35 | 1 | 35 | 35 | 35 | 35 |
| C10H22 | 75 | 1 | 75 | 75 | 75 | 75 |
| regular_d3_n6 | 2 | 1 | 1 | 2 | 2 | 2 |
| regular_d3_n8 | 5 | 1 | 1 | 5 | 5 | 5 |
| regular_d3_n10 | 18 | 1 | 1 | 18 | 18 | 18 |
| regular_d3_n12 | 78 | 1 | 1 | 78 | 78 | 78 |
| regular_d4_n8 | 6 | 1 | 1 | 6 | 6 | 6 |
| regular_d4_n10 | 54 | 1 | 1 | 54 | 54 | 54 |
| strong_regular_S16 | 2 | 1 | 1 | 1 | 2 | 2 |

## Atom-order invariance

Each molecule was renumbered five times; a non-zero count means the ID depends on input atom order.

| algorithm | molecules | variants | molecules with false splits |
| --- | ---: | ---: | ---: |
| formula | 38 | 228 | 0 |
| wl1 | 38 | 228 | 0 |
| wl2 | 38 | 228 | 0 |
| canonical_smiles | 38 | 228 | 0 |
| inchi | 38 | 228 | 0 |

## Timing scaling

| family | method | log-log slope |
| --- | --- | ---: |
| chain | wl1_full | 1.894 |
| chain | wl1_k4 | 0.954 |
| chain | canonical_smiles | 1.49 |
| chain | inchi | 0.991 |
| regular3 | wl1_full | 0.739 |
| regular3 | wl1_k4 | 0.73 |
| regular3 | canonical_smiles | 1.184 |
| regular3 | inchi | 0.66 |
| regular3 | wl2_full | 3.19 |

## Brute-force audit (n <= 8 heavy atoms)

| group | molecules | brute force classes | canonical SMILES classes | VF2 classes | agree |
| --- | ---: | ---: | ---: | ---: | --- |
| C2H6O | 2 | 2 | 2 | 2 | True |
| C3H6O | 6 | 6 | 6 | 6 | True |
| C3H8O | 3 | 3 | 3 | 3 | True |
| C4H10O | 7 | 7 | 7 | 7 | True |
| C4H8 | 5 | 5 | 5 | 5 | True |
| C4H6 | 7 | 7 | 7 | 7 | True |
| C5H12 | 3 | 3 | 3 | 3 | True |
| C6H6 | 5 | 5 | 5 | 5 | True |
| C5H5NO | 2 | 2 | 2 | 2 | True |
| regular_d3_n6 | 2 | 2 | 2 | 2 | True |
| regular_d3_n8 | 5 | 5 | 5 | 5 | True |

## Notes

Standard InChI merges 2-pyridone and 2-hydroxypyridine (same InChI). Use non-isomeric canonical SMILES, or fixed-H InChI, when constitutional identity must be exact.

1-WL is exact on trees and on all non-regular graphs in these sets, but it is provably incomplete: every pair of non-isomorphic d-regular graphs with the same n and d has the same 1-WL hash. 2-WL fixes the 6-vertex K(3,3)/prism case but still collides on the Shrikhande and 4x4 rook graphs (6-regular S16). Only canonical labelling is exact.
