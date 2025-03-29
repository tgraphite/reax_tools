# ReaxTools 1.1

(This doucment is translated by Claude-3.7 from Chinese version, be careful)

A high-performance C++ tool for analyzing molecular dynamics trajectories, with a focus on reactive systems.

[中文版说明](README_zh.md)

## Updates in 1.1
- Added reaction flow graph functionality to track reactions between molecules
- Improved molecule identification algorithm for better accuracy
- Support for both .xyz and .lammpstrj file formats
- Fixed multiple stability issues

## Features

- **High Performance**: Completely rewritten in C++, optimized for speed and memory efficiency
- **Trajectory Analysis**: Direct processing of trajectory files (.xyz, .lammpstrj)
- **Species Identification**: Identifies molecular species using van der Waals radii
- **Reaction Tracking**: Tracks reactions between frames and builds reaction flow graphs
- **Visualization**: Generates DOT format graphs for reaction visualization
- **Low Dependencies**: Only requires the C++ standard library and fmt

## Usage

```
Usage: 
  -f .xyz/.lammpstrj file -> [TRAJ MODE] determine molecules by van der Waals radius
  -s lammps reaxff/species file (species.out) -> [SPECIES MODE] determine species by file input

[TRAJ MODE SETTINGS]
  -r radius scaling factor (default 1.2)
  -t type names, split in comma, e.g. C,H,O,N,S,F

[SPECIES MODE SETTINGS]
  -me merge molecules into groups by element number (default C), i.e. mols have 1~4 Carbons -> group_C1-C4
  -mr merge range for the process above, split in comma (default: 1,4,8,16)
  -rc rescale group weight by selected atom number, not mol number, i.e. C2H4 -> weight 2, C4H8 -> weight 4 (default: no)
  --order output formulas in correct element order, split in comma (default: C,H,O,N,S,F,P)
```

### Examples

Process an XYZ trajectory with carbon, hydrogen, oxygen, and nitrogen atoms:
```bash
./reax_tools -f trajectory.xyz -t C,H,O,N
```

Process a LAMMPS trajectory with custom van der Waals radius scaling:
```bash
./reax_tools -f trajectory.lammpstrj -t C,H,O,N -r 1.3
```

## Output Files

- `.species.csv` - Contains molecular species statistics and evolution data
- `.dot` - Reaction flow graph that can be visualized with tools like Graphviz

## Test Output
```
Frame: 1 Atoms: 12902 Bonds: 10559 Mols: 5173
Frame: 2 Atoms: 12902 Bonds: 10223 Mols: 5180
Frame: 3 Atoms: 12902 Bonds: 10181 Mols: 5243
Frame: 4 Atoms: 12902 Bonds: 10073 Mols: 5284
Frame: 5 Atoms: 12902 Bonds: 10031 Mols: 5326
C1 : 86 74 95 108 104 
C102H8O44 : 0 0 0 1 0 
C104H8O43 : 0 0 1 0 0 
C10H1O3 : 0 0 0 0 1 
C10H1O4 : 0 0 1 0 1 
C10H1O5 : 0 0 0 0 1 
(...)

Save file ../examples/FeCHON_5frames.csv successfully.
=== Reaction Flow Report ===
Total nodes (species): 40
Total edges (reactions): 49

Top reactions:
1: C2H1 -> C1H1 (count: 17)
2: C1H1 -> C2H1 (count: 13)
3: C3H2 -> C2H1 (count: 9)
4: C2H1 -> C3H2 (count: 6)
5: C2O1 -> C1O1 (count: 5)
6: C2O1 -> C2 (count: 4)
7: C1H1O1 -> C2H2O1 (count: 4)
8: H2N1 -> H3N1 (count: 4)
9: C2H2 -> C3H3 (count: 3)
10: C3H3 -> C2H2 (count: 3)
Reaction flow graph saved to ../examples/FeCHON_5frames.dot
```

## Implementation Details

- Tick-Tock alternating frame reading to save memory and support frame-to-frame analysis
- K-D tree algorithm for efficient neighbor atom searching
- Chemical bond determination based on van der Waals radii
- Depth-first search algorithm for building molecular graphs
- Molecule similarity calculation based on atom ID intersection/union to track reactions

## Performance

- Test case: 350MB lammpstrj trajectory, 13,000 atoms, 1,000 frames
- Performance: ~160 seconds on a single-core processor, less than 100MB memory usage
- Species file mode is even faster, typically taking only a few seconds

## Future Development Plans

- Reaction chains: Track the complete lifecycle of molecules from creation to destruction
- Group transfer: Analyze the migration of molecular fragments between different molecules
- Reaction kinetics: Calculate reaction rates and reaction orders
- Visualization: Provide richer graphical output and reporting functionality

## Building from Source

```bash
mkdir build
cd build
cmake ..
make -j4
```

## License

No License, do anything you want.