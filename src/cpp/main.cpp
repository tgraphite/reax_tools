#include <chrono>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "argparser.h"
#include "fmt/core.h"
#include "reax_counter.h"
#include "smiles.h"
#include "string_tools.h"
#include "universe.h"

static void write_run_log(double elapsed_sec,
                          const std::string& command_line,
                          const Universe& uv) {
    FILE* fp = create_file("reax_tools.log");
    if (fp == nullptr) return;

    fmt::print(fp, "[run]\n");
    fmt::print(fp, "schema_version = 2.1\n");
    fmt::print(fp, "program = ReaxTools\n");
    fmt::print(fp, "version = {}\n", "2.3");
    fmt::print(fp, "identity_model = canonical-smiles-v1\n");
    fmt::print(fp, "smiles_model = reax-canonical-v1\n");
    fmt::print(fp, "command = {}\n", command_line);
    fmt::print(fp, "input = {}\n", INPUT_FILE);
    fmt::print(fp, "output_dir = {}\n", OUTPUT_DIR);
    fmt::print(fp, "elapsed_s = {:.6f}\n", elapsed_sec);
    fmt::print(fp, "frames = {}\n", uv.frames_processed);
    fmt::print(fp, "atoms_first_frame = {}\n", uv.atoms_first_frame);
    fmt::print(fp, "atoms_last_frame = {}\n\n", uv.atoms_last_frame);

    fmt::print(fp, "[parameters]\n");
    fmt::print(fp, "vdw_scale = {:.6f}\n", RVDW_FACTOR);
    fmt::print(fp, "threads = {}\n", NUM_THREADS);
    fmt::print(fp, "max_neighbors = {}\n", MAX_NEIGH);
    fmt::print(fp, "rings = {}\n", FLAG_NO_RINGS ? "disabled" : "enabled");
    fmt::print(fp, "reaction_events = {}\n", FLAG_TRACK_REACTIONS ? "enabled" : "disabled");
    fmt::print(fp, "transfer_flow = {}\n", FLAG_NO_REACTIONS ? "disabled" : "enabled");
    fmt::print(fp, "stable_time = {}\n", STABLE_TIME_FRAMES);
    fmt::print(fp, "timestep_fs = {:.6f}\n", TIMESTEP_FS);
    fmt::print(fp, "sampling_freq = {}\n", SAMPLING_FREQ);
    fmt::print(fp, "smiles_max_atoms = {}\n\n", reax_smiles::kMaxSmilesAtoms);

    fmt::print(fp, "[results]\n");
    fmt::print(fp, "species = {}\n", uv.species_counter->species_count());
    fmt::print(fp, "formulas = {}\n", uv.species_counter->formula_count());
    fmt::print(fp, "events = {}\n", uv.reaction_tracker ? uv.reaction_tracker->get_event_count() : 0);
    fmt::print(fp, "transfer_edges = {}\n\n",
               uv.reaction_tracker ? uv.reaction_tracker->get_transfer_edge_count() : 0);

    fmt::print(fp, "[outputs]\n");
    fmt::print(fp, "species_count.csv\n");
    fmt::print(fp, "bond_count.csv\n");
    fmt::print(fp, "atom_bonded_num_count.csv\n");
    fmt::print(fp, "ring_count.csv\n");
    if (FLAG_TRACK_REACTIONS) {
        fmt::print(fp, "reaction_events.csv\n");
        fmt::print(fp, "reaction_event_pairs.csv\n");
    }
    if (!FLAG_NO_REACTIONS) fmt::print(fp, "transfer_flow.csv\n");
    fmt::print(fp, "molecules.csv\n");
    fmt::print(fp, "reax_tools.log\n\n");

    fmt::print(fp, "[validation]\n");
    fmt::print(fp, "audit_basis = raw atom-overlap transfer events\n");
    fmt::print(fp, "status = generated audit-ready raw outputs\n");
    fclose(fp);
}

int run_analysis(int argc, const char* const* argv) {
    auto start_time = std::chrono::high_resolution_clock::now();
    fmt::print("ReaxTools 2.3\n");
    fmt::print("High-performance reactive MD post-processing\n\n");

    std::string command_line;
    for (int i = 0; i < argc; ++i) {
        if (i > 0) command_line += " ";
        command_line += argv[i];
    }

    ArgParser parser = init_argparser();
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help" || arg == "help") {
            parser.print_help();
            return 0;
        }
    }

    try {
        parser.parse_args(argc, argv);
        parser.operate_for_all();
    }
    catch (const std::exception& e) {
        std::cerr << "ArgParser error: " << e.what() << std::endl;
        return 1;
    }

    if (!std::filesystem::exists(OUTPUT_DIR)) {
        std::filesystem::create_directory(OUTPUT_DIR);
    }

    Universe uv;
    uv.process_traj();

    uv.bond_counter->save_file("bond_count.csv");
    uv.ring_counter->save_file("ring_count.csv");
    uv.atom_bonded_num_counter->save_file("atom_bonded_num_count.csv");

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed_sec = std::chrono::duration<double>(end_time - start_time).count();
    write_run_log(elapsed_sec, command_line, uv);

    fmt::print("\nSummary\n");
    fmt::print("  elapsed time     {:.3f} s\n", elapsed_sec);
    fmt::print("  output directory {}\n", OUTPUT_DIR);
    fmt::print("  run log          reax_tools.log\n");

    return 0;
}

int main(int argc, char** argv) {
    return run_analysis(argc, argv);
}
