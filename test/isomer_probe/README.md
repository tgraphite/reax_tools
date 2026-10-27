# Isomer identity probe (ReaxTools 2.3)

This probe checks structural identity on a pair of constitutional isomers with
the same molecular formula, C2H6O:

- ethanol: C0-C1-O2
- dimethyl ether (DME): C0-O2-C1

The atom order (C, C, O, H x 6) and the atom IDs are identical in both
frames, so the trajectory is a valid, atom-conserving isomerization.

## Files

| File | Contents |
| --- | --- |
| `isomerization_ethanol_to_dme.xyz` | 2 frames, 9 atoms: ethanol then DME |
| `mixture_ethanol_plus_dme.xyz` | 1 frame, 18 atoms: both isomers coexist |
| `out_isomerization/` | 2.3 analysis of the 2-frame trajectory |
| `out_mixture/` | 2.3 analysis of the coexisting-isomer frame |

The geometries were generated once with RDKit/ETKDG+MMFF and are stored as
plain XYZ; RDKit is not needed to reproduce the analysis.

## Reproduce

```bash
./bin/reax_tools_core -f test/isomer_probe/isomerization_ethanol_to_dme.xyz \
  -o test/isomer_probe/out_isomerization -nt 1
./bin/reax_tools_core -f test/isomer_probe/mixture_ethanol_plus_dme.xyz \
  -o test/isomer_probe/out_mixture -nt 1
```

## Observed result (ReaxTools 2.3)

Both frames contain 9 atoms, 8 bonds and 1 molecule. The two structures receive
different IDs and different canonical SMILES:

```csv
id,formula,smiles
9888429903061026892,C2H6O,[C]([H])([H])([C]([H])([H])([H]))[O]([H])
12567411285834319415,C2H6O,[C]([H])([H])([H])[O][C]([H])([H])([H])
```

The isomerization is now a 1 -> 1 reaction event (same formula, different
structure ID):

```csv
event_id,frame,n_reactants,n_products,reactant_ids,product_ids,reactant_formulas,product_formulas,reactant_smiles,product_smiles
0,2,1,1,9888429903061026892,12567411285834319415,C2H6O,C2H6O,...
```

When both isomers coexist in one frame, each is counted separately:

```csv
frame,9888429903061026892,12567411285834319415
1,1,1
```

`species_count.csv` is a wide table: the first column is `frame`, and every
other column is a structure ID. The formula and SMILES are in `molecules.csv`.

## Interpretation

The SMILES are connectivity-only (all bonds single) and every H is explicit,
so a parser cannot add implicit hydrogens. The identity separates constitutional
isomers; it intentionally does not encode bond orders, charge, spin, or
stereochemistry. Components above 200 atoms fall back to formula-level identity.
