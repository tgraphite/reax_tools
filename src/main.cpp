#include <iostream>

#include "defines.h"
#include "getopt.h"
#include "reax_species.h"
#include "string_tools.h"
#include "universe.h"

int main(int argc, char *argv[]) {
    std::string mode;
    std::string traj_file;
    std::string species_file;
    std::string xyz_file;
    std::string merge_target = "C";

    std::string usage =
        "Usage: \n\
	-f <.xyz/.lammpstrj file> analyze trajectory, build speceis and reaction flows.\n\
	-s <lammps reaxff/species file (spieces.out)> clean this file only.\n\n\
    \
	[TRAJ ANALYSIS SETTINGS]\n\
    -t <type_names>, mandatory when using lammpstrj file, split in comma, e.g. C,H,O,N,S,F\n\
	-r <radius_scaling_factor> (default 1.2)\n\
	-nt <num_threads> (default 4)\n\
	--dump write lammps data for each frame (may cost more time)\n\
    \
	[SPECIES ANALYSIS SETTINGS]\n\
    -me if merge molecules into groups by element number (default C), i.e. mols have 1~4 Carbons -> group_C1-C4\n\
	-mr <merge_ranges>, split in comma (default: 1,4,8,16)\n\
	-rc rescale group weight by selected atom number, not mol number, i.e. C2H4 -> weight 2, C4H8 -> weight 4 (default: no)\n\
	--order output formulas in correct element order, split in comma (default: C,H,O,N,S,F,P)\n\
    \
    [EXAMPLES]\n\
    reax_tools -f traj.lammpstrj -t C,H,O,N,S,F -r 1.2 -nt 4 -me C -rc --dump\n\
    reax_tools -s species.out -me C -mr 1,4,8,16 -rc\n\
    ";

    std::vector<std::string> type_names;
    std::vector<int> merge_range = {1, 4, 8, 16};
    std::vector<std::string> sort_order = sorted_elements;

    int num_threads = 4;
    float rvdw_scale = 1.2f;
    bool if_merge_by_element = false;
    bool if_merge_rescale = false;
    bool if_dump_lammps_data = false;
    bool if_dump_bond_count = true;
    int max_reactions = 50;

    std::map<std::string, int> opts_nvals = {{"-f", 1},  {"-s", 1},     {"-r", 1},      {"-t", 1}, {"-me", 1},
                                             {"-mr", 1}, {"-rc", 1},    {"--order", 1}, {"-h", 0}, {"--help", 0},
                                             {"-nt", 1}, {"--dump", 0}, {"--topo", 1}};

    std::map<std::string, std::vector<std::string>> opts_vals = neo_getopt(argc, argv, opts_nvals, usage);

    for (auto &pair : opts_vals) {
        std::string opt = pair.first;
        std::vector<std::string> vals = pair.second;

        if (opt == "-f") {
            traj_file = vals[0];
            mode = "traj";
        } else if (opt == "-s") {
            species_file = vals[0];
            mode = "species";
        } else if (opt == "--topo") {
            mode = "topo";
            xyz_file = vals[0];
        } else if (opt == "-r") {
            rvdw_scale = std::stof(vals[0]);
        } else if (opt == "-t") {
            type_names = split(vals[0], ",");
        } else if (opt == "-nt") {
            num_threads = std::stoi(vals[0]);
        } else if (opt == "-me") {
            if_merge_by_element = true;
            merge_target = vals[0];
        } else if (opt == "-mr") {
            std::vector<std::string> merge_range_str = split(vals[0], ",");
            merge_range.clear();
            for (auto &val : merge_range_str) {
                merge_range.push_back(std::stoi(val));
            }
        } else if (opt == "-rc") {
            if_merge_rescale = true;
        } else if (opt == "--order") {
            sort_order.clear();
            sort_order = split(vals[0], ",");
        } else if ((opt == "-h") || (opt == "--help")) {
            std::cout << usage << std::endl;
            exit(0);
        } else if (opt == "--dump") {
            if_dump_lammps_data = true;
        }
    }

    // Process, read lammpstrj mode
    if (mode == "traj") {
        Universe uv;

        uv.process_traj(traj_file, type_names, rvdw_scale, num_threads, if_dump_lammps_data);

        if (if_merge_by_element) {
            uv.reax_species->merge_by_element(merge_target, merge_range, if_merge_rescale);
        }

        uv.reax_species->rename_all_formulas(sort_order);
        uv.reax_species->brief_report();
        uv.reax_species->save_file(traj_file);

        uv.reax_flow->brief_report();
        uv.reax_flow->save_graph(traj_file, max_reactions);
    }

    // Read lammps species.out mode
    else if (mode == "species") {
        ReaxSpecies reax_species(species_file);
        if (if_merge_by_element) {
            reax_species.merge_by_element(merge_target, merge_range, if_merge_rescale);
        }
        reax_species.rename_all_formulas(sort_order);
        reax_species.brief_report();
        reax_species.save_file();
    }

    // Create topo file mode
    // else if (mode == "topo") {
    //     System sys;
    //     sys.load_file(xyz_file);
    //     sys.search_neigh(2.5 * rvdw_scale, 10);
    //     sys.build_molecules();
    //     for (auto &mol : sys.mols) {
    //         std::cout << mol->info() << std::endl;
    //     }
    // }
}