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
    std::string merge_target = "C";

    std::string usage =
        "Usage: \n\
	-f .xyz/.lammpstrj file -> [TRAJ MODE] determine molecules by van der Waals radius\n\
    -nt num_threads (default 4)\n\
	-s lammps reaxff/species file (spieces.out) -> [SPECIES MODE] determine species by file input\n\n\
	[TRAJ MODE SETTINGS]\n\
	-r raidus scaling factor (default 1.2)\n\
	-t type names, split in comma, e.g. C,H,O,N,S,F\n\
	-b begin frame (default 1)\n\
	-e end frame (default: all frames)\n\
	--step frame step (default 1)\n\n\
	[SPECIES MODE SETTINGS]\n\
	-me merge molecules into groups by element number (default C), i.e. mols have 1~4 Carbons -> group_C1-C4\n\
	-mr merge range for the process above, split in comma (default: 1,4,8,16)\n\
	-rc rescale group weight by selected atom number, not mol number, i.e. C2H4 -> weight 2, C4H8 -> weight 4 (default: no)\n\
	--order output formulas in correct element order, split in comma (default: C,H,O,N,S,F,P)\n";

    std::vector<std::string> type_names;
    std::vector<int> merge_range = {1, 4, 8, 16};
    std::vector<std::string> sort_order = sorted_elements;

    int num_threads = 4;
    float rvdw_scale = 1.2f;
    bool if_merge_by_element = false;
    bool if_merge_rescale = false;

    std::map<std::string, int> opts_nvals = {{"-f", 1},  {"-s", 1},     {"-r", 1},  {"-t", 1},
                                             {"-me", 1}, {"-mr", 1},    {"-rc", 1}, {"--order", 1},
                                             {"-h", 0},  {"--help", 0}, {"-nt", 1}};

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
        } else if ((opt == "-h") or (opt == "--help")) {
            std::cout << usage << std::endl;
            exit(0);
        }
    }

    // Process, read lammpstrj mode
    if (mode == "traj") {
        Universe uv;

        uv.process_traj(traj_file, type_names, rvdw_scale, num_threads);

        if (if_merge_by_element) {
            uv.reax_species->merge_by_element(merge_target, merge_range, if_merge_rescale);
        }

        uv.reax_species->rename_all_formulas(sort_order);
        uv.reax_species->brief_report();
        uv.reax_species->save_file(traj_file);

        uv.reax_flow->brief_report();
        uv.reax_flow->save_graph(traj_file);

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
}