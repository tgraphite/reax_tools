# energetic_v4 (2.1) vs energetic_v5 (2.2)

Input and parameters are the same: `test/energetic.xyz`, `vdw_scale = 1.2`,
40 frames, reaction tracking enabled. v5 was produced by the 2.2 core; v4 is
the historical formula-level output.

Reproduce:

```bash
python3 analysis/compare_v4_v5.py --v4 test/energetic_v4 --v5 test/energetic_v5
```

## Bottom line

- **Manual merge of v5 by formula reproduces v4 species counts exactly.**
  126 formulas x 40 frames: **0 mismatched cells**.
- `bond_count.csv`, `atom_bonded_num_count.csv` and `ring_count.csv` are
  byte-identical.
- The v5 formula set is identical to v4; v5 simply resolves **183 structures**
  instead of 126 formulas.
- Reaction events: v5 is a **strict superset** of v4 at formula level.
  **0 missing**, **2 extra** events. Both extras are the same-formula
  tautomerization of H2N8O11 at frames 36 and 37.
- After old-style common-species cancellation, the net reaction signature set
  is identical; the only difference is +1 `<cancelled>` at frames 36 and 37.
- Transfer flow: 984 v4 formula edges vs 985 v5 collapsed edges. **All counts
  and atom-transfer sums match except one new self-loop**:
  `H2N8O11 -> H2N8O11`, count 2, atom transfer 42. This is the collapsed
  footprint of the two tautomerization events.

## Counts

| item | v4 | v5 | v5 merged by formula |
| --- | ---: | ---: | ---: |
| species rows / structures | 126 formulas | 183 structures | 126 formulas |
| species-count mismatches | - | - | 0 |
| isomers per formula | 1 for all | 86 x 1, 29 x 2, 7 x 3, 3 x 4, 1 x 6 | 1 |
| raw reaction events | 1415 | 1417 | 1417 |
| unique raw event signatures | 1241 | 1243 | 1241 |
| after common-species cancellation | 1221 | 1221 | identical key set; +2 `<cancelled>` counts |
| transfer edges | 984 | 985 | 984 + 1 self-loop |

## The only substantive difference: H2N8O11

The formula H2N8O11 has two structures in v5:

- A (ID 8330235418410714537): frames 15-35 and 37-40
- B (ID 17611810749496458707): frame 36 only

The two 1 -> 1 events are A -> B at frame 36 and B -> A at frame 37. They have
the same 21 atoms (atom_overlap = 21 in `reaction_event_pairs.csv`).

The heavy-atom skeleton is the same; the structures differ by **one hydrogen
position**: A has an `NH2` group, B has `NH` plus `OH`. This is a
tautomerization/proton transfer, not a heavy-atom rearrangement. The old
formula-level trigger could not see it because the formula does not change.

When the two v5 flow edges (A->B and B->A) are collapsed to formulas, they
produce the single formula self-loop `H2N8O11 -> H2N8O11` (count 2,
atom_transfer 42). v4 has no such edge.

## Interpretation

- For **composition/counting**, v5 + manual formula merge is exactly v4.
- For **net formula reactions** after cancelling common species, v5 is also
  exactly v4; the two extras cancel to no net reaction.
- For **raw event and transfer tables**, v5 intentionally contains two more
  events and one more self-loop. These are real same-formula structural
  changes that the 2.1 identity model was blind to, not a regression.
- If strict v4-compatible output is required, drop 1->1 events whose reactant
  and product formulas are identical (or apply the old common-species
  cancellation) before comparing; then the event signature set is identical.
  The formula self-loop must be dropped alongside them.
