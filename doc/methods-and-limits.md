# Analysis semantics and current limits

Source review: 2.2 development. This describes the current implementation, not an accuracy certification.

## From coordinates to events

1. `System::load_xyz/load_lammpstrj` reads coordinates and, where present, orthorhombic box lengths. Four-column XYZ assigns IDs by row order. The reader has special-case layouts, not a general extended-XYZ `Properties` parser. LAMMPS timesteps are read past rather than retained as the output time axis.
2. `System::process_this` constructs pair cutoffs `r_ij = scale * (R_i + R_j)/2`. Defaults come from `argparser.cpp`; some element values are generic placeholders. Periodic neighbor search uses cell lists; non-periodic search uses pair comparisons. The periodic neighbor candidate cap is 10 by default. Dense environments need particular care because this cap can truncate candidates before distance sorting.
3. `System::build_bonds_by_radius` sorts each atom's candidates by normalized squared distance, greedily accepting edges while both coordination caps permit. This is not a global graph optimization. Atom ordering/ties can affect ambiguous environments. Connected components define molecules. A subsequent heuristic increases bond orders using residual valence; these are not electronic bond orders.
4. `Molecule::update_topology` computes a canonical SMILES from the perceived element/H graph (all bonds treated as single; no charge or spin perception) and hashes it with 64-bit FNV-1a. Components above 200 atoms skip SMILES and use a 64-bit hash of the formula, so their SMILES field is empty and isomers are not separated. A hash collision is possible in principle and is reported when the same ID is seen with different SMILES.
5. `ReactionTracker::process_frame` detects atoms whose owning molecule's **structure ID** changed. A closure over the old and new molecules groups those seeds into events. This resembles traversal of the overlap graph between consecutive molecule partitions; the comment “symmetric difference” should not be interpreted as a proof of a distinct new graph algorithm. A pure bond-order change is invisible because SMILES is connectivity-only.
6. Each event records positive intersections between reactant and product atom-ID sets. If `a(e,r,p)` is that intersection size, a species edge has `count = sum 1[a>0]` and `atom_transfer = sum a` over matching event/instance pairs. The event may be many-to-many; a projected edge is not a balanced elementary reaction. Same-formula pairs inside a larger detected event may create self-loops.
7. Python filters and plots the raw data. The default network displays selected species/transfers. Event plots can cancel common formulas; that is a presentation transformation, not temporal recrossing suppression. `flow` additionally suppresses some cyclic/lateral information.

## What the audits establish

`tools/check_reax_outputs.py` checks expected files, selected headers, species-ID references and optional historical counts. `tools/audit_event_network_consistency.py` independently aggregates the pair table and compares it to the flow table. The latter does not read the event index, raw coordinates or atom-ID membership, so it cannot by itself detect omitted reactions, orphan event references or chemically wrong bonds. Both files can agree while sharing an upstream error.

`ReactionEvent::check_atom_conservation` compares the union of atom IDs on both sides and warns on mismatch; it does not reject the event. Pair files store intersection sizes, not the complete intersecting atom-ID lists. Coordinates and bond-change lists are not written for events.

## Limits to record with a calculation

| Topic | Current behavior / requirement |
| --- | --- |
| Periodicity | Orthorhombic, all axes periodic if a cell is present. Off-diagonal lattice terms and per-axis periodicity flags are not applied. A small cell shorter than the neighbor cutoff also needs validation. |
| Identity | Canonical SMILES (64-bit ID) for components of at most 200 atoms; constitutional isomers are separated. SMILES is connectivity-only: bond orders, charge, spin, stereochemistry and tautomer normalization are not represented. Components above 200 atoms fall back to formula-level identity. |
| Time filtering | `--stable-time` is parsed; `should_mark_stable` is defined but not called by the raw-event path. No HMM or lifetime-based raw-event suppression is implemented. |
| Time axis | CSV frame IDs start at 1. Convert to physical time from the actual saved-frame interval and starting time; do not assume the internal defaults recover original MD timestamps. |
| Ring counts | Bounded graph-ring search, defaults 5–8 atoms. Not an ice classifier, aromaticity assignment, or validated universal primitive-ring census. |
| Memory | Frames are batched, but tracked molecules, events, counters and wide species tables accumulate. Total memory is not independent of trajectory length. |
| Reproducibility | Event IDs/order and representative structures can vary with traversal/threading. Compare normalized event multisets including multiplicity, and retain each run's matching pair IDs. |
| Metadata | Run manifest says `canonical-smiles-v1` / `reax-canonical-v1`. Custom radius/valence maps and all timing options are not fully captured. Save the full command and input checksum externally. |
| Empty/disabled results | Event/pair/flow files are written with headers even when no events are detected; disabling reaction tracking omits them. Always use a new output directory to avoid stale files from an earlier run. |
| Scope of a network | Counts and repeated atom-overlap sums over observed events. No rate constants, activation barriers, elementary mechanism, or guaranteed complete reaction catalogue. |

## Reproducibility record

Keep the source commit, compiler/build flags, full command, input checksum, format and atom-ID convention, cell/periodicity, MD engine and potential, saved-frame time interval, radius and coordination overrides, thread count, output checksums, and all plotting filters. For a methods paper, also record an independently labelled reaction set, sensitivity results and timing/memory protocols. See [paper preparation](paper-preparation/README.md).

Implementation anchors: [`system.cpp`](../src/cpp/system.cpp), [`system.h`](../src/cpp/system.h), [`smiles.cpp`](../src/cpp/smiles.cpp), [`cell_list.cpp`](../src/cpp/cell_list.cpp), [`reax_track.cpp`](../src/cpp/reax_track.cpp), [`universe.cpp`](../src/cpp/universe.cpp), [`main.cpp`](../src/cpp/main.cpp).
