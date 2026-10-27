# Isomer identity probe

This probe checks the species-identity model of ReaxTools 2.1 on a pair of
constitutional isomers with the same molecular formula, `C2H6O`:

- ethanol: `C0-C1-O2` (8 bonds, 1 molecule)
- dimethyl ether (DME): `C0-O2-C1` (8 bonds, 1 molecule)

The atom order (`C, C, O, H x 6`) and the atom IDs are identical in both
frames, so the trajectory is a valid, atom-conserving isomerization.

## Files

| File | Contents |
| --- | --- |
| `isomerization_ethanol_to_dme.xyz` | 2 frames, 9 atoms: ethanol then DME |
| `mixture_ethanol_plus_dme.xyz` | 1 frame, 18 atoms: both isomers coexist |
| `out_isomerization/` | analysis of the 2-frame trajectory |
| `out_mixture/` | analysis of the coexisting-isomer frame |

The geometries were generated once with RDKit/ETKDG+MMFF and are stored as
plain XYZ; RDKit is not needed to reproduce the analysis.

## Reproduce

```bash
./bin/reax_tools_core -f test/isomer_probe/isomerization_ethanol_to_dme.xyz \
  -o test/isomer_probe/out_isomerization -nt 1
./bin/reax_tools_core -f test/isomer_probe/mixture_ethanol_plus_dme.xyz \
  -o test/isomer_probe/out_mixture -nt 1
```

## Observed result (ReaxTools 2.1, commit dc6237b)

Both frames contain 9 atoms, 8 bonds and 1 molecule, and the connectivity
changes from ethanol to DME, but the identity is unchanged:

```
frame,molecule_id,formula,count      # out_isomerization/species_count.csv
1,1124119999,C2H6O,1
2,1124119999,C2H6O,1
```

No `reaction_events.csv` is written: the 1 -> 1 isomerization is invisible to
the event trigger because `Molecule::hash` is `FNV-1a(formula)`.

When both isomers coexist in one frame, they are counted as two copies of one
species:

```
frame,molecule_id,formula,count      # out_mixture/species_count.csv
1,1124119999,C2H6O,2
```

`molecules.json` contains a single `C2H6O` node whose `example` graph is
ethanol; the DME graph is not retained.

## Interpretation

This is expected behaviour for the documented `formula-hash-v2` identity
model, not a regression. It is not a structural identifier: constitutional
isomers, conformers, stereoisomers and charge/spin states with the same
formula all share one ID. See the "Raw Outputs" section of the main README.
