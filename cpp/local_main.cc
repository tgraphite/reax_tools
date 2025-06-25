#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "argparser.h"
#include "defines.h"
#include "fmt/core.h"
#include "reax_species.h"
#include "string_tools.h"
#include "universe.h"

int main(int argc, char **argv) {
    ArgParser parser("reax_tools_local", "ReaxFF Trajectory Analyzer (Local Test)");
    parser.add_argument("--traj", "-f", "Analyze trajectory file (.xyz/.lammpstrj)", "file", "", true, false,
                        "path/to/file");
    parser.add_argument("--types", "-t", "element types splitted in commas", "Traj analysis", "", false, false,
                        "str,str,...(e.g. C,H,O,N...)");
    parser.add_argument("--radius", "-r", "scaling factor of vdW radii", "Traj analysis", "1.2", false, false, "float");
    parser.add_argument("--threads", "-nt", "number of threads", "Traj analysis", "4", false, false, "int");
    parser.add_argument("--dump", "", "dump lammps data file (.data) for each frame", "Traj analysis", "", false, true);
    parser.add_argument("--reaxflow-threshold", "--reaxflow", "similarity threshold for analyze reaction flows",
                        "Traj analysis", "0.25", false, false, "float");
    parser.add_argument("--reduce-reactions", "-rr", "reduce reverse reactions", "Traj analysis", "false", false, true);
    parser.add_argument("--merge-element", "-me", "merge species groups by an element type", "Species analysis", "C",
                        false, false, "e.g. C");
    parser.add_argument("--merge-ranges", "-mr", "merge group range, split by commas", "Species analysis", "1,4,8,16",
                        false, false, "e.g. 1,4,8,16");
    parser.add_argument("--rescale-count", "-rc",
                        "rescale group weight by actual numbers of atoms (instead of molecules)", "Species analysis",
                        "", false, true);
    parser.add_argument("--element-order", "--order", "set element order of outputting formulas", "Species analysis",
                        "", false, false, "e.g. C,H,O,N,S,F,P");

    if (!parser.parse_args(argc, argv)) {
        return 1;
    }

    std::string traj_file = parser.has_option("--traj") ? parser.get<std::string>("--traj") : "";

    std::cout << "traj_file: [" << traj_file << "]" << std::endl;
    if (!std::filesystem::exists(traj_file)) {
        std::cerr << "Error: File does not exist: " << traj_file << std::endl;
        return 1;
    }

    if (traj_file.empty()) {
        std::cerr << "Error: Must define input file (-f or --traj)" << std::endl;
        return 1;
    }

    std::vector<std::string> type_names;
    if (parser.has_option("--types")) {
        type_names = parser.get<std::vector<std::string>>("--types");
    }
    float rvdw_scale = parser.get<float>("--radius");
    int num_threads = parser.get<int>("--threads");
    float reaxflow_threshold = parser.get<float>("--reaxflow-threshold");
    bool if_dump_lammps_data = parser.has_flag("--dump");
    bool if_reduce_reactions = parser.has_flag("--reduce-reactions") || parser.has_flag("-rr");
    std::string merge_target = parser.has_option("--merge-element") ? parser.get<std::string>("--merge-element") : "";
    std::vector<int> merge_range = parser.get<std::vector<int>>("--merge-ranges");
    bool if_merge_by_element = !merge_target.empty();
    bool if_merge_rescale = parser.has_flag("--rescale-count") || parser.has_flag("-rc");
    std::vector<std::string> sort_order;
    if (parser.has_option("--element-order")) {
        sort_order = parser.get<std::vector<std::string>>("--element-order");
    } else {
        sort_order = default_order;
    }
    int max_reactions = 60;

    if (ends_with(traj_file, "lammpstrj") && type_names.empty()) {
        std::cerr << "Error: Must define element types when using lammpstrj file. (--types)" << std::endl;
        return 1;
    }

    Universe uv;
    std::string output_dir = "output/";

    fmt::print("=== Reax Tools Trajectory Analysis (Local) ===\n");
    fmt::print("Input file: {}, output dir: {}\n", traj_file, output_dir);

    if (!std::filesystem::exists(output_dir)) {
        std::filesystem::create_directory(output_dir);
    } else {
        for (const auto &entry : std::filesystem::directory_iterator(output_dir)) {
            std::filesystem::remove_all(entry);
        }
    }

    if (std::find(type_names.begin(), type_names.end(), "X") != type_names.end()) {
        fmt::print("**** NOTE **** : You are using the X flag to ignore this type of atoms.\n");
        fmt::print("**** NOTE **** : This may cause inaccuracy when X is not an elementary substance.\n");
    }

    uv.process_traj(traj_file, output_dir, type_names, rvdw_scale, num_threads, if_dump_lammps_data,
                    reaxflow_threshold);

    if (if_merge_by_element) {
        uv.reax_species->merge_by_element(merge_target, merge_range, if_merge_rescale);
    }

    uv.reax_species->brief_report();
    uv.reax_species->save_file_to_dir(output_dir);
    uv.reax_flow->brief_report();
    uv.reax_flow->save_graph(output_dir, max_reactions, true, if_reduce_reactions);
    uv.reax_flow->dump_smiles(output_dir);
    uv.reax_flow->draw_molecules(output_dir);

    return 0;
}
