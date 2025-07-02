#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "defines.h"
#include "fmt/core.h"
#include "reax_species.h"
#include "string_tools.h"
#include "universe.h"

Universe::Universe() {}

Universe::~Universe() {
    if (system != nullptr) {
        delete system;
        system = nullptr;
    }

    if (last_system != nullptr) {
        delete last_system;
        last_system = nullptr;
    }

    if (reax_species != nullptr) {
        delete reax_species;
        reax_species = nullptr;
    }

    if (reax_flow != nullptr) {
        delete reax_flow;
        reax_flow = nullptr;
    }
}

void Universe::flush() {
    // only for fallbacked serial.
    if (last_system != nullptr) {
        delete last_system;
        last_system = nullptr;
    }
    last_system = std::move(system);
    system = nullptr;
}

/**
 * Fallback to serial in wasm. num_threads is ignored.
 */
void Universe::process_traj(std::string& file_path, std::string& output_dir, std::vector<std::string>& type_names,
                            const float& rvdw_scale, const int& num_threads, const bool& if_dump_lammps_data,
                            const bool& if_no_reax_flow) {
    int curr_frame_id = 1;
    int max_neigh = 10;
    bool is_first_frame = true;

    std::ifstream file(file_path);
    std::string bond_count_filepath = output_dir + "bond_count.csv";
    std::string ring_count_filepath = output_dir + "ring_count.csv";
    std::string atom_bonded_num_count_filepath = output_dir + "atom_bonded_num_count.csv";

    // The highest calling stack, only do this once.
    while (file.is_open() && !file.eof()) {
        if (curr_frame_id > 1) {
            flush();
        }

        system = new System();

        system->set_frame_id(curr_frame_id);
        system->set_types(type_names);
        system->set_max_neigh(max_neigh);
        system->set_rvdw_scale(rvdw_scale);
        system->set_reax_flow(this->reax_flow);
        system->set_reax_species(this->reax_species);

        if (ends_with(file_path, ".lammpstrj"))
            system->load_lammpstrj(file);
        else if (ends_with(file_path, ".xyz"))
            system->load_xyz(file);
        else {
            std::cerr << "Using unsupported format, check your file and suffix!";
            exit(1);
        }

        if (system->atoms.size() == 0) {
            continue;
        };
        curr_frame_id++;

        system->set_prev_sys(last_system);
        system->process_this();
        system->process_reax_species();
        system->process_reax_flow();

        if (if_dump_lammps_data) {
            std::string lammps_data_file = output_dir + "frame_" + std::to_string(system->frame_id) + ".data";
            system->dump_lammps_data(lammps_data_file);
        }

        if (system->frame_id == 1) {
            fmt::print("Atom Types: ");
            for (auto& pair : system->type_itos) {
                fmt::print("{}: {}, ", pair.first, pair.second);
            }
            fmt::print("\n");
            fmt::print("\n");

            is_first_frame = true;
        } else {
            is_first_frame = false;
        }

        system->dump_bond_count(bond_count_filepath, is_first_frame);
        system->dump_ring_count(ring_count_filepath, is_first_frame);
        system->dump_atom_bonded_num_count(atom_bonded_num_count_filepath, is_first_frame);
        system->finish();
    }
    fmt::print("\n\n");
    reax_species->analyze_frame_formulas();
}
