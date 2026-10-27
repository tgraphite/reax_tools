#include <chrono>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "argparser.h"
#include "fmt/core.h"
#include "reax_counter.h"
#include "string_tools.h"
#include "universe.h"

static void write_run_log(double elapsed_sec) {
    FILE* fp = create_file("reax_tools.log");
    if (fp == nullptr) return;

    fmt::print(fp, "[run]\n");
    fmt::print(fp, "program = ReaxTools\n");
    fmt::print(fp, "version = {}\n", "2.2");
    fmt::print(fp, "identity_model = canonical-smiles-v1\n");
    fmt::print(fp, "smiles_model = reax-canonical-v1\n");
    fmt::print(fp, "input = {}\n", INPUT_FILE);
    fmt::print(fp, "output_dir = {}\n", OUTPUT_DIR);
    fmt::print(fp, "elapsed_s = {:.6f}\n\n", elapsed_sec);

    fmt::print(fp, "[parameters]\n");
    fmt::print(fp, "vdw_scale = {:.6f}\n", RVDW_FACTOR);
    fmt::print(fp, "threads = {}\n", NUM_THREADS);
    fmt::print(fp, "max_neighbors = {}\n", MAX_NEIGH);
    fmt::print(fp, "rings = {}\n", FLAG_NO_RINGS ? "disabled" : "enabled");
    fmt::print(fp, "reaction_events = {}\n", FLAG_TRACK_REACTIONS ? "enabled" : "disabled");
    fmt::print(fp, "transfer_flow = {}\n\n", FLAG_NO_REACTIONS ? "disabled" : "enabled");

    fmt::print(fp, "[outputs]\n");
    fmt::print(fp, "species_count.csv\n");
    fmt::print(fp, "bond_count.csv\n");
    fmt::print(fp, "atom_bonded_num_count.csv\n");
    fmt::print(fp, "ring_count.csv\n");
    fmt::print(fp, "reaction_events.csv\n");
    fmt::print(fp, "reaction_event_pairs.csv\n");
    fmt::print(fp, "transfer_flow.csv\n");
    fmt::print(fp, "molecules.csv\n");
    fmt::print(fp, "reax_tools_manifest.json\n");
    fmt::print(fp, "reax_tools.log\n\n");

    fmt::print(fp, "[validation]\n");
    fmt::print(fp, "audit_basis = raw atom-overlap transfer events\n");
    fmt::print(fp, "status = generated audit-ready raw outputs\n");
    fclose(fp);
}

static void write_manifest(double elapsed_sec) {
    FILE* fp = create_file("reax_tools_manifest.json");
    if (fp == nullptr) return;

    fmt::print(fp, "{{\n");
    fmt::print(fp, "  \"schema_version\": \"2.0\",\n");
    fmt::print(fp, "  \"program\": \"ReaxTools\",\n");
    fmt::print(fp, "  \"version\": \"2.2\",\n");
    fmt::print(fp, "  \"identity_model\": \"canonical-smiles-v1\",\n");
    fmt::print(fp, "  \"smiles_model\": \"reax-canonical-v1\",\n");
    fmt::print(fp, "  \"input\": {{\"trajectory\": \"{}\"}},\n", INPUT_FILE);
    fmt::print(fp, "  \"output_directory\": \"{}\",\n", OUTPUT_DIR);
    fmt::print(fp, "  \"elapsed_s\": {:.6f},\n", elapsed_sec);
    fmt::print(fp, "  \"parameters\": {{\n");
    fmt::print(fp, "    \"vdw_scale\": {:.6f},\n", RVDW_FACTOR);
    fmt::print(fp, "    \"threads\": {},\n", NUM_THREADS);
    fmt::print(fp, "    \"ring_detection\": {},\n", FLAG_NO_RINGS ? "false" : "true");
    fmt::print(fp, "    \"reaction_events\": {},\n", FLAG_TRACK_REACTIONS ? "true" : "false");
    fmt::print(fp, "    \"transfer_flow\": {}\n", FLAG_NO_REACTIONS ? "false" : "true");
    fmt::print(fp, "  }},\n");
    fmt::print(fp, "  \"files\": [\n");
    fmt::print(fp, "    \"species_count.csv\",\n");
    fmt::print(fp, "    \"bond_count.csv\",\n");
    fmt::print(fp, "    \"atom_bonded_num_count.csv\",\n");
    fmt::print(fp, "    \"ring_count.csv\",\n");
    fmt::print(fp, "    \"reaction_events.csv\",\n");
    fmt::print(fp, "    \"reaction_event_pairs.csv\",\n");
    fmt::print(fp, "    \"transfer_flow.csv\",\n");
    fmt::print(fp, "    \"molecules.csv\",\n");
    fmt::print(fp, "    \"reax_tools.log\",\n");
    fmt::print(fp, "    \"reax_tools_manifest.json\"\n");
    fmt::print(fp, "  ]\n");
    fmt::print(fp, "}}\n");
    fclose(fp);
}

int run_analysis(int argc, const char* const* argv) {
    auto start_time = std::chrono::high_resolution_clock::now();
    fmt::print("ReaxTools 2.2\n");
    fmt::print("High-performance reactive MD post-processing\n\n");

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
    write_run_log(elapsed_sec);
    write_manifest(elapsed_sec);

    fmt::print("\nSummary\n");
    fmt::print("  elapsed time     {:.3f} s\n", elapsed_sec);
    fmt::print("  output directory {}\n", OUTPUT_DIR);
    fmt::print("  manifest         reax_tools_manifest.json\n");

    return 0;
}

int main(int argc, char** argv) {
    return run_analysis(argc, argv);
}
