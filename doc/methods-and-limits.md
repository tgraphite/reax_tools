# Analysis semantics and current limits

Source review: 2026-09-11, commit `dc6237b4a67a5e3a2f11322c78b3fe96193b7aef`. This describes the current implementation, not an accuracy certification.

## From coordinates to events

1. `System::load_xyz/load_lammpstrj` reads coordinates and, where present, orthorhombic box lengths. Four-column XYZ assigns IDs by row order. The reader has special-case layouts, not a general extended-XYZ `Properties` parser. LAMMPS timesteps are read past rather than retained as the output time axis.
2. `System::process_this` constructs pair cutoffs `r_ij = scale * (R_i + R_j)/2`. Defaults come from `argparser.cpp`; some element values are generic placeholders. Periodic neighbor search uses cell lists; non-periodic search uses pair comparisons. The periodic neighbor candidate cap is 10 by default. Dense environments need particular care because this cap can truncate candidates before distance sorting.
3. `System::build_bonds_by_radius` sorts each atom's candidates by normalized squared distance, greedily accepting edges while both coordination caps permit. This is not a global graph optimization. Atom ordering/ties can affect ambiguous environments. Connected components define molecules. A subsequent heuristic increases bond orders using residual valence; these are not electronic bond orders.
4. `Molecule::update_topology` hashes the formula using 32-bit FNV-1a. Formula-identical components share species identity regardless of connectivity, charge or spin. Hash collisions are also possible in principle. A representative graph in `molecules.json` must not be treated as the unique structure of that species.
5. `ReactionTracker::process_frame` detects atoms whose owning molecule's **formula hash** changed. A closure over the old and new molecules groups those seeds into events. This resembles traversal of the overlap graph between consecutive molecule partitions; the comment “symmetric difference” should not be interpreted as a proof of a distinct new graph algorithm. Same-formula atom exchange with no changed-formula seed is invisible to this trigger.
6. Each event records positive intersections between reactant and product atom-ID sets. If `a(e,r,p)` is that intersection size, a species edge has `count = sum 1[a>0]` and `atom_transfer = sum a` over matching event/instance pairs. The event may be many-to-many; a projected edge is not a balanced elementary reaction. Same-formula pairs inside a larger detected event may create self-loops.
7. Python filters and plots the raw data. The default network displays selected species/transfers. Event plots can cancel common formulas; that is a presentation transformation, not temporal recrossing suppression. `flow` additionally suppresses some cyclic/lateral information.

## What the audits establish

`tools/check_reax_outputs.py` checks expected files, selected headers, species-ID references and optional historical counts. `tools/audit_event_network_consistency.py` independently aggregates the pair table and compares it to the flow table. The latter does not read the event index, raw coordinates or atom-ID membership, so it cannot by itself detect omitted reactions, orphan event references or chemically wrong bonds. Both files can agree while sharing an upstream error.

`ReactionEvent::check_atom_conservation` compares the union of atom IDs on both sides and warns on mismatch; it does not reject the event. Pair files store intersection sizes, not the complete intersecting atom-ID lists. Snapshot packages provide more detailed atom data for a subset of event arities.

## Limits to record with a calculation

| Topic | Current behavior / requirement |
| --- | --- |
| Periodicity | Orthorhombic, all axes periodic if a cell is present. Off-diagonal lattice terms and per-axis periodicity flags are not applied. A small cell shorter than the neighbor cutoff also needs validation. |
| Identity | Formula level; no isomer, charge, spin or stereochemical resolution. Same-formula exchanges can be omitted. |
| Time filtering | `--stable-time` is parsed; `should_mark_stable` is defined but not called by the raw-event path. No HMM or lifetime-based raw-event suppression is implemented. |
| Time axis | CSV frame IDs start at 1. Convert to physical time from the actual saved-frame interval and starting time; do not assume the internal defaults recover original MD timestamps. |
| Snapshot time | Tracked molecules copy coordinates when created; snapshots use those stored instances. Reactant geometry may predate the event boundary. Validate against adjacent original frames before using a snapshot as a transition illustration. |
| Ring counts | Bounded graph-ring search, defaults 5–8 atoms. Not an ice classifier, aromaticity assignment, or validated universal primitive-ring census. |
| Memory | Frames are batched, but tracked molecules, events, counters and snapshots accumulate. Total memory is not independent of trajectory length. |
| Reproducibility | Event IDs/order and representative structures can vary with traversal/threading. Compare normalized event multisets including multiplicity, and retain each run's matching pair IDs. |
| Metadata | Run manifest says `formula-hash-v2`; molecule metadata still says `formula-hash-v1`. Custom radius/valence maps and all timing options are not fully captured. Save the full command and input checksum externally. |
| Empty/disabled results | Some files are omitted when events are empty or features disabled although the manifest lists them. Always use a new output directory to avoid stale files from an earlier run. |
| Scope of a network | Counts and repeated atom-overlap sums over observed events. No rate constants, activation barriers, elementary mechanism, or guaranteed complete reaction catalogue. |

## Reproducibility record

Keep the source commit, compiler/build flags, full command, input checksum, format and atom-ID convention, cell/periodicity, MD engine and potential, saved-frame time interval, radius and coordination overrides, thread count, output checksums, and all plotting filters. For a methods paper, also record an independently labelled reaction set, sensitivity results and timing/memory protocols. See [paper preparation](paper-preparation/README.md).

Implementation anchors: [`system.cpp`](../src/cpp/system.cpp), [`system.h`](../src/cpp/system.h), [`cell_list.cpp`](../src/cpp/cell_list.cpp), [`reax_track.cpp`](../src/cpp/reax_track.cpp), [`universe.cpp`](../src/cpp/universe.cpp), [`main.cpp`](../src/cpp/main.cpp).
