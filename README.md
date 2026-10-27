# ReaxTools 2.1

ReaxTools is a local post-processing toolkit for reactive molecular dynamics trajectories from LAMMPS, GPUMD, CP2K, and related workflows.

Version 2.1 is a major cleanup release. The C++ core now writes stable raw audit files, while Python owns filtering, plotting, and presentation. The atom-transfer table is auditable from raw reaction events; `plot` and `network` provide compact molecule-transfer and reaction-transfer diagrams.

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
- The Python plotting installation includes PyGraphviz and PyVista/VTK. Network plots need Graphviz; building PyGraphviz may also require Graphviz development headers and Python development headers (on Debian/Ubuntu: `graphviz libgraphviz-dev python3-dev`).

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
| `species_count.csv`         | Species count table: `frame,molecule_id,formula,count`               |
| `bond_count.csv`            | Bond-type counts by frame                                            |
| `atom_bonded_num_count.csv` | Atom coordination counts by frame                                    |
| `ring_count.csv`            | Ring counts by frame                                                 |
| `reaction_events.csv`       | Lightweight reaction-event index: frame, arity, hashes, and formulas |
| `reaction_event_pairs.csv`  | Per-event atom-overlap pairs used to audit atom transfer             |
| `transfer_flow.csv`         | Aggregated atom transfer between molecular species                   |
| `reaction_snapshots/`       | C++ reaction snapshot packages for montage rendering                 |
| `molecules.json`            | Molecule identity records and example graph definitions              |
| `reax_tools.log`            | Run metadata and audit notes                                         |
| `reax_tools_manifest.json`  | Structured file manifest for Python and web clients                  |

Species IDs are deterministic 32-bit FNV-1a hashes of the molecular formula. Isomers share an ID; this is not a structural identifier. `molecules.json` contains representative graphs, not an exhaustive isomer catalogue. Current run manifests say `formula-hash-v2`, while the molecule writer still labels its metadata `formula-hash-v1`; historical example outputs also use the older label/IDs. Keep outputs from each run together.

The table describes the default analysis with detected events. Zero-event runs can omit event, pair, transfer, and snapshot files; the manifest currently lists anticipated filenames even when absent. The audit commands below target populated default runs and will flag missing files. Disabling reaction analysis also changes available outputs.

## Audited Transfer

The reaction-transfer analysis is based on atom ownership transfer. `reaction_events.csv` is a compact event index, while `reaction_event_pairs.csv` records per-event reactant-product atom overlaps. `transfer_flow.csv` is exactly the aggregation of `reaction_event_pairs.csv`. For single- and bimolecular snapshot rendering, the C++ core writes self-contained `.rxtsnap.json` packages under `reaction_snapshots/`.

Run the built-in audit:

```bash
python3 tools/check_reax_outputs.py test/energetic_v4 \
  --expect-events 1415 \
  --expect-transfer-edges 984 \
  --expect-self-loop-count 115

python3 tools/audit_event_network_consistency.py test/energetic_v4 \
  --expect-event-pairs 4695
```

`energetic_v4` has the event headers expected by the current schema checker. The older `energetic_v3` still passes pair/flow aggregation, but its event headers are incompatible with that checker. Do not overwrite the historical fixtures when running your own analysis.

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
reax_tools plot -f test/energetic_v3
```

Individual products:

```bash
reax_tools counts -f test/energetic_v3
reax_tools counts -f test/energetic_v3 -t H2O,N2,H3N
reax_tools events -f test/energetic_v3 --max 15
reax_tools network -f test/energetic_v3 --max-reactions 60
reax_tools snapshots -f test/energetic_v3
reax_tools flow -f test/energetic_v3
reax_tools focus -f test/energetic_v3 --centers H3N N3O4 H2O
```

`plot` and the default `network` command write two transfer diagrams:

- `molecule_transfer_top24.png`: the 24 most active molecular species and their direct transfers.
- `reaction_transfer_top48.png`: the 48 strongest reaction transfers and the molecular species involved.

`network` also accepts explicit chemical filters:

- `--max-reactions N`: keep the strongest N reaction transfers.
- `--max-molecules N`: keep the strongest N molecular species and the direct transfers among them.
- `--max-subgraphs N`: keep the largest connected molecule groups.

These filters combine with AND semantics. `flow` remains an explicit experimental command for one-way role-layer summaries.

## Energetic Example

The bundled energetic fixture comes from an ammonium-dinitramide ReaxFF run at 373 K. Its 40 stored frames span 0.2–78.2 ps at 2 ps intervals, as checked against the source trajectory and run log. It is a software regression example, not a validated high-temperature decomposition mechanism. The plots show formula counts and atom transfers involving `H4N`, `N3O4`, `H2O`, and `N2`; geometric bond perception and sampling affect those observations. The images below are retained historical outputs in `test/energetic_v3`, not a new benchmark. See [provenance](doc/paper-preparation/local-cases.md).

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
