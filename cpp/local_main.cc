#include <chrono>
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

int main(int argc, char** argv) {

    auto start_time = std::chrono::high_resolution_clock::now();

    ArgParser parser("reax_tools", "Reactive MD Trajectory Analyzer (Local machine)");

    //  parser.add_argument(const std::string& name,
    //                      const std::string& short_name,
    //                      const std::string& description,
    //                      const std::string& group = "General",
    //                      const std::string& default_value = "",
    //                      bool required = false,
    //                      bool is_flag = false,
    //                      const std::string& value_type = "",
    //                      std::function<bool(const std::string&)> validator = nullptr)

    parser.add_argument("--traj", "-f", "Analyze trajectory file (.xyz/.lammpstrj)", "Input file", "", true, false,
                        "path/to/file");
    parser.add_argument("--types", "-t", "element types splitted in commas", "Key input parameters", "", false, false,
                        "str,str,...(e.g. C,H,O,N...)");
    parser.add_argument("--radius", "-r", "scaling factor of vdW radii", "Key input parameters", "1.2", false, false,
                        "float");
    parser.add_argument("--type-radius", "-tr", "set atomic radius for element, e.g. -tr N:1.5", "Key input parameters",
                        "", false, false, "Element:Radius");

    parser.add_argument("--threads", "-nt", "number of threads", "Performance", "4", false, false, "int");

    parser.add_argument("--dump", "", "dump lammps data file (.data) for each frame", "Misc", "", false, true);
    parser.add_argument("--no-reactions", "", "disable reaction analysis", "Network output options", "", false, true);
    parser.add_argument("--no-reduce-reactions", "-norr", "disable reduce reverse reactions", "Network output options",
                        "false", false, true);
    parser.add_argument("--max-reactions", "", "set max reactions in default network", "Network output options", "60",
                        false, false);

    parser.add_argument("--merge-element", "-me", "merge species groups by an element type", "Species output options",
                        "C", false, false, "e.g. C");
    parser.add_argument("--merge-ranges", "-mr", "merge group range, split by commas", "Species output options",
                        "1,4,8,16", false, false, "e.g. 1,4,8,16");
    parser.add_argument("--rescale-count", "-rc",
                        "rescale group weight by actual numbers of atoms (instead of molecules)",
                        "Species output options", "", false, true);
    parser.add_argument("--element-order", "--order", "set element order of outputting formulas",
                        "Species output options", "", false, false, "e.g. C,H,O,N,S,F,P");

    if (!parser.parse_args(argc, argv)) {
        std::cerr << "Error: Failed to parse arguments. Check your input." << std::endl;
        return 1;
    }

    std::string traj_file = parser.get<std::string>("--traj");
    if (!std::filesystem::exists(traj_file)) {
        std::cerr << "Error: File does not exist: " << traj_file << std::endl;
        return 1;
    }

    std::vector<std::string> type_names;
    if (parser.has_option("--types")) {
        type_names = parser.get<std::vector<std::string>>("--types");
    }

    int num_threads = parser.get<int>("--threads");
    float rvdw_scale = parser.get<float>("--radius");

    bool if_dump_lammps_data = parser.has_flag("--dump");
    bool if_no_reax_flow = parser.has_flag("--no-reactions");
    bool if_no_reduce_reactions = parser.has_flag("--no-reduce-reactions") || parser.has_flag("-norr");
    int max_reactions = parser.get<int>("--max-reactions");

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
        for (const auto& entry : std::filesystem::directory_iterator(output_dir)) {
            std::filesystem::remove_all(entry);
        }
    }

    if (std::find(type_names.begin(), type_names.end(), "X") != type_names.end()) {
        fmt::print("**** NOTE **** : You are using the X flag to ignore this type of atoms.\n");
        fmt::print("**** NOTE **** : This may cause inaccuracy when X is not an elementary substance.\n");
    }

    if (parser.has_option("--type-radius")) {
        std::vector<std::string> type_radius_list = parser.get<std::vector<std::string>>("--type-radius");
        for (const auto& tr : type_radius_list) {
            auto pos = tr.find(':');
            if (pos != std::string::npos) {
                std::string elem = tr.substr(0, pos);
                float radius = std::stof(tr.substr(pos + 1));
                default_atomic_radius[elem] = radius; // update global map
            } else {
                std::cerr << "Invalid --type-radius format: " << tr << std::endl;
            }
        }
    }

    uv.process_traj(traj_file, output_dir, type_names, rvdw_scale, num_threads, if_dump_lammps_data, if_no_reax_flow);

    if (if_merge_by_element) {
        uv.reax_species->merge_by_element(merge_target, merge_range, if_merge_rescale);
    }

    uv.reax_species->brief_report();
    uv.reax_species->save_file_to_dir(output_dir);
    uv.reax_flow->brief_report();
    uv.reax_flow->save_graph(output_dir, max_reactions, true, if_no_reduce_reactions);
    uv.reax_flow->dump_smiles(output_dir);
    uv.reax_flow->draw_molecules(output_dir);

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed_sec = std::chrono::duration<double>(end_time - start_time).count();
    fmt::print("Total elapsed time: {:.3f} seconds\n", elapsed_sec);

    return 0;
}
