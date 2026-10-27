# ReaxTools 2.3

ReaxTools is a local post-processing toolkit for reactive molecular dynamics trajectories from LAMMPS, GPUMD, CP2K, and related workflows.

Version 2.3 keeps **canonical-SMILES structural identity** and removes the JSON manifest: `reax_tools.log` is the single run record. `plot` also writes top-48 2D structure montages (`key_molecules_page*.png`) using abundance and reaction-graph degree. The C++ core writes raw audit files, while Python owns labels, filtering, plotting, and presentation.

## Install

```bash
git clone https://github.com/tgraphite/reax_tools
cd reax_tools
bash install_reax_tools.sh
```

The installer builds `bin/reax_tools_core`, creates a project-local Python virtual environment (`venv/`), installs the plotting dependencies from `bin/requirements.txt`, and adds `bin/` to your shell startup file. Restart the terminal or source that file after installation. The `reax_tools` command automatically runs inside the venv when it exists. If venv creation fails, the installer attempts a user-level pip installation.

Requirements:

- CMake 3.14+ and a C++17 compiler (GCC 9+ recommended; GCC 8.x is also supported and links `libstdc++fs` automatically).
- The installer accepts Python 3.8+ and searches for `python3.13` down to `python3.8`, then `python3`. The unpinned plotting dependencies have their own Python requirements; this is not a tested compatibility matrix for every accepted Python version.
- The Python plotting installation includes PyGraphviz. Network plots need Graphviz; building PyGraphviz may also require Graphviz development headers and Python development headers (on Debian/Ubuntu: `graphviz libgraphviz-dev python3-dev`).
- RDKit is **optional**: install it only if you want `reax_tools molecules` to draw structures from SMILES (`pip install rdkit`).

The C++ core does not require RDKit, Boost, or Graphviz. The Python plotting layer has additional dependencies listed above.

You can also build the C++ core directly without the Python layer:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./bin/reax_tools_core -f test/energetic.xyz -o reax_tools_output
```

## Quick Start

Analyze a trajectory:

```bash
reax_tools -f test/energetic.xyz -o reax_tools_output
```

`reax_tools -f ...` is a shortcut for analyze mode. The explicit form is also supported:

```bash
reax_tools analyze -f trajectory.xyz -o reax_tools_output
```

Generate the default plots from an output directory:

```bash
reax_tools plot -f reax_tools_output
```

Inspect available commands:

```bash
reax_tools
reax_tools analyze --help
reax_tools network --help
```

## Input Notes

ReaxTools supports:

- XYZ with `element x y z` atom rows; selected extended XYZ layouts with a `Lattice="..."` comment are supported. The reader does not implement arbitrary `Properties` column layouts. Preserve atom order across frames when IDs are implicit.
- Plain XYZ is treated as non-periodic when no cell is supplied. CP2K coordinate-only XYZ therefore needs cell information added before periodic analysis.
- LAMMPS dump files with the `.lammpstrj` suffix, stable `id`, numeric `type`, and either `x y z` or `xs ys zs` columns. Use `-t C,H,O,N` to map types 1, 2, 3, 4 to elements.

Periodic analysis currently uses orthorhombic box lengths and periodicity along all three axes. Triclinic tilt and per-axis boundary flags are not implemented. Converting a skew cell to XYZ alone does not resolve this limitation. Keep the atom population and IDs fixed across frames. See [method and input limits](doc/methods-and-limits.md) before analyzing a new trajectory format.

Useful analyze options:

```bash
reax_tools -f traj.xyz -o out
reax_tools -f dump.lammpstrj -t C,H,O,N -o out
reax_tools -f traj.xyz -r 1.1
reax_tools -f traj.xyz -tr N:1.5 O:1.6
reax_tools -f traj.xyz -tv N:4 O:3
reax_tools -f traj.xyz --no-rings
```

`-r` sets the distance cutoff scale (default 1.2): a pair is eligible when `distance < 0.5 * (radius_i + radius_j) * scale`. Candidate neighbors are then selected subject to element-specific coordination caps. `-tr` overrides radii and `-tv` overrides caps. These are geometric heuristics, not ReaxFF bond orders; validate them for your chemistry. Increasing the cutoff admits more candidates but does not guarantee monotonic changes after coordination limiting.

`--stable-time` is currently a legacy parameter: it does not filter the raw reaction events. Raw events can contain short-lived bond fluctuations. Frame numbers start at 1 and refer to stored trajectory frames, not the original MD step numbers.

## Raw Outputs

The C++ core writes raw, audit-oriented files:

| File                        | Meaning                                                              |
| --------------------------- | -------------------------------------------------------------------- |
| `species_count.csv`         | Wide count table: `frame,<id>,<id>,...` (one column per structure)   |
| `bond_count.csv`            | Bond-type counts by frame                                            |
| `atom_bonded_num_count.csv` | Atom coordination counts by frame                                    |
| `ring_count.csv`            | Ring counts by frame                                                 |
| `reaction_events.csv`       | Event index with IDs, formulas, and SMILES for reactants/products    |
| `reaction_event_pairs.csv`  | Per-event atom-overlap pairs used to audit atom transfer             |
| `transfer_flow.csv`         | Aggregated atom transfer between structures, with formulas and SMILES|
| `molecules.csv`             | Identity table: `id,formula,smiles`                                 |
| `key_molecules_page*.png`   | Optional 2D structure montages for the top-48 important structures  |
| `reax_tools.log`            | Single human-readable run metadata and audit record                  |

Species IDs are deterministic **64-bit FNV-1a hashes of canonical SMILES**. A different connectivity gives a different SMILES and therefore a different ID, even when the formula is identical. Formal charges are emitted only when needed to keep standard valence parseable (N valence 4, O valence 3); they are not a charge-equilibration model. Molecules with more than 200 atoms are treated at formula level: no SMILES is generated and the ID is a hash of the formula; `molecules.csv` leaves the SMILES field empty. This cutoff keeps kerogen/polymer-scale components and accidental whole-system molecules tractable; within the cutoff, constitutional isomers are separated.

The core writes the event/pair/flow CSV headers even when no events are detected (with reaction tracking enabled), so downstream readers see a stable file set. Disabling reaction analysis omits those files.

## Audited Transfer

The reaction-transfer analysis is based on atom ownership transfer. `reaction_events.csv` is a compact event index, while `reaction_event_pairs.csv` records per-event reactant-product atom overlaps. `transfer_flow.csv` is exactly the aggregation of `reaction_event_pairs.csv`. Reaction events now carry IDs, formulas, and SMILES; no coordinates or bond-change lists are written.

Run the built-in audit:

```bash
./bin/reax_tools_core -f test/energetic.xyz -o /tmp/reax_check
python3 tools/check_reax_outputs.py /tmp/reax_check
python3 tools/audit_event_network_consistency.py /tmp/reax_check
```

The checked-in `test/energetic_v3` and `test/energetic_v4` directories are historical 2.1 outputs with the old formula-level schema. They are not compatible with the 2.3 checkers; regenerate with the current core before auditing.

For the locally available `gpumd_v3` fixture, the pair/flow aggregation audit also passes (its source trajectory and some supporting files are not tracked in Git):

```text
reaction events: 6522
reaction event pairs: 72471
transfer edges: 18391
atom transfer total: 1733312
```

These are historical fixture checks, not expected counts for every input or version. The pair aggregation audit checks equality of edge counts and atom-overlap sums; it does not independently validate the chemistry or prove that every physical reaction was detected. In particular, same-formula rearrangements and exchanges can be missed by the current event trigger.

An edge's `count` is the number of overlapping reactant-product instance pairs, and `atom_transfer` sums their shared atoms. Neither is a rate constant, mass flux, or count of distinct atoms over the whole trajectory. A raw multi-reactant event can contribute several edges. Self-loops are retained in raw data; plotting filters may remove them.

## Plotting Commands

Default plots:

```bash
reax_tools plot -f /tmp/reax_check
```

Individual products:

```bash
reax_tools counts -f /tmp/reax_check
reax_tools counts -f /tmp/reax_check -t H2O,N2,H3N
reax_tools events -f /tmp/reax_check --max 15
reax_tools network -f /tmp/reax_check --max-reactions 60
reax_tools flow -f /tmp/reax_check
reax_tools focus -f /tmp/reax_check --centers H3N N3O4 H2O
# Requires optional RDKit:
reax_tools molecules -f /tmp/reax_check --species H2O
```

`plot` also writes `key_molecules_page1.png` through `key_molecules_page3.png`: the 48 most important structures by time-averaged abundance and reaction-graph degree, drawn as 2D structures with `formula (a)` labels. Drawings apply a **display-only** residual-valence greedy bond-order assignment (the same idea as the core bond builder); the stored SMILES and identity remain all-single and connectivity-only. This requires RDKit; if RDKit is absent, `plot` skips the montage with a note. Species count pagination is unchanged.

`plot` and the default `network` command write two transfer diagrams:

- `molecule_transfer_top24.png`: the 24 most active molecular species and their direct transfers.
- `reaction_transfer_top48.png`: the 48 strongest reaction transfers and the molecular species involved.

`network` also accepts explicit chemical filters:

- `--max-reactions N`: keep the strongest N reaction transfers.
- `--max-molecules N`: keep the strongest N molecular species and the direct transfers among them.
- `--max-subgraphs N`: keep the largest connected molecule groups.

These filters combine with AND semantics. `flow` remains an explicit experimental command for one-way role-layer summaries.

## Energetic Example

The bundled energetic fixture comes from an ammonium-dinitramide ReaxFF run at 373 K. Its 40 stored frames span 0.2–78.2 ps at 2 ps intervals, as checked against the source trajectory and run log. It is a software regression example, not a validated high-temperature decomposition mechanism. The plots show species counts and atom transfers involving `H4N`, `N3O4`, `H2O`, and `N2`; geometric bond perception and sampling affect those observations. The images below are retained historical 2.1 outputs in `test/energetic_v3` (formula-level identity), not 2.3 outputs or a new benchmark. See [provenance](doc/paper-preparation/local-cases.md).

Species counts:

![Species count](test/energetic_v3/species_count.png)

Reaction-event summary:

![Reaction events](test/energetic_v3/reaction_events.png)

Molecule-transfer view:

![Molecule transfer](test/energetic_v3/molecule_transfer_top24.png)

Reaction-transfer view:

![Reaction transfer](test/energetic_v3/reaction_transfer_top48.png)

## Plot Style Templates

Python plotting reads YAML settings in this order:

1. Command-line `--template`.
2. `reax_tools_plot.yaml` or `reax_tools_template.yaml` in the output directory.
3. Bundled `src/python/reax_tools_viz/default_plot_template.yaml`.

This lets the C++ raw output remain stable while plot style evolves in Python.

## Repository Layout

```text
reax_tools/
├── install_reax_tools.sh       # one-step local install
├── bin/
│   ├── reax_tools              # user-facing command router
│   ├── reax_tools_core         # C++ analysis core
│   └── requirements.txt        # Python plotting dependencies
├── src/
│   ├── cpp/                    # raw analysis core
│   └── python/reax_tools_viz/  # plotting and filtering layer
├── tools/                      # audit utilities
└── test/                       # small fixtures and v2.1 example outputs
```

## Citation

The software does not yet have a methods paper. Cite the repository and record the version/commit and access date. [Methods and limitations](doc/methods-and-limits.md) describe the current analysis semantics. Internal [paper preparation](doc/paper-preparation/README.md) records the evidence, candidate cases, and proposed outline; the [paper-freeze archive](doc/paper-freeze/README.md) records the reproducibility metadata and pending publication gates. Neither is a manuscript.

If ReaxTools helps your work, please cite or acknowledge:

```text
Hanxiang Chen. ReaxTools: A high performance ReaxFF/AIMD/MLP-MD post-process code [Computer software].
https://github.com/tgraphite/reax_tools
```
