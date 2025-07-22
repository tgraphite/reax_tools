#include <chrono>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "argparser.h"
#include "fmt/core.h"
#include "reax_counter.h"
#include "string_tools.h"
#include "universe.h"

#ifndef WASM_MODE
int main(int argc, char** argv)
#else
extern "C" int cpp_main(int argc, const char** argv)
#endif
{
    auto start_time = std::chrono::high_resolution_clock::now();
    fmt::print("ReaxTools: High performance reactive MD post-processer\n");

    ArgParser parser = init_argparser();
    try {
        parser.parse_args(argc, argv);
        parser.operate_for_all();
    } catch (const std::exception& e) {
        std::cerr << "ArgParser error: " << e.what() << std::endl;
        return 1;
    }

    std::string output_dir = "reax_tools_output/";
    if (!std::filesystem::exists(output_dir)) {
        std::filesystem::create_directory(output_dir);
    }

    Universe uv;
    uv.process_traj();

    uv.reax_counter->brief_report();
    uv.reax_counter->save_file();

    uv.reax_flow->brief_report();
    uv.reax_flow->save_graph();

#ifndef WASM_MODE
    uv.reax_flow->dump_smiles();
    uv.reax_flow->draw_molecules();
#endif

    uv.bond_counter->save_file("bond_count.csv");
    uv.ring_counter->save_file("ring_count.csv");
    uv.atom_bonded_num_counter->save_file("atom_bonded_num_count.csv");

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed_sec = std::chrono::duration<double>(end_time - start_time).count();
    fmt::print("Total elapsed time: {:.3f} seconds\n", elapsed_sec);

    return 0;
}
