#include "universe.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#include "atom.h"
#include "defines.h"
#include "fmt/format.h"
#include "molecule.h"
#include "reax_species.h"
#include "string_tools.h"

Universe::Universe() {
    last_system = nullptr;
    current_systems.clear();
}

Universe::~Universe() {
    for (auto sys : current_systems) {
        if (sys != nullptr) {
            delete sys;
        }
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
    if (last_system != nullptr) {
        delete last_system;
        last_system = nullptr;
    }

    if (current_systems.size() > 0 && current_systems[current_systems.size() - 1] != nullptr) {
        last_system = current_systems[current_systems.size() - 1];
        current_systems[current_systems.size() - 1] = nullptr;
    } else {
        last_system = nullptr;
    }

    for (auto &sys : current_systems) {
        if (sys != nullptr) {
            delete sys;
        }
    }
    current_systems.clear();
}

void Universe::process_traj(std::string &file_path, std::string &output_dir, std::vector<std::string> &type_names,
                            const float &rvdw_scale, const int &num_threads, const bool &if_dump_lammps_data,
                            const float &reaxflow_threshold) {
    int curr_frame_id = 1;
    int max_neigh = 10;
    float neigh_radius = 2.8 * rvdw_scale;
    bool is_first_frame = true;

    std::ifstream file(file_path);
    std::string bond_count_filepath = output_dir + "bond_count.csv";
    std::string ring_count_filepath = output_dir + "ring_count.csv";

    // The highest calling stack, only do this once.

    while (file.is_open() && !file.eof()) {
        if (curr_frame_id > 1) {
            flush();
        }

        // Serial read. systems readed in the right order.
        for (int thread_id = 0; thread_id < num_threads; thread_id++) {
            System *curr_sys = new System();

            curr_sys->set_frame_id(curr_frame_id);
            curr_sys->set_types(type_names);
            curr_sys->set_max_neigh(max_neigh);
            curr_sys->set_rvdw_scale(rvdw_scale);
            curr_sys->set_reax_flow(this->reax_flow);
            curr_sys->set_reax_species(this->reax_species);
            curr_sys->set_reaxflow_threshold(reaxflow_threshold);

            current_systems.push_back(curr_sys);

            if (ends_with(file_path, ".lammpstrj"))
                curr_sys->load_lammpstrj(file);
            else if (ends_with(file_path, ".xyz"))
                curr_sys->load_xyz(file);

            if (curr_sys->atoms.size() == 0) {
                continue;
            };
            curr_frame_id++;
        }

        for (int thread_id = 0; thread_id < num_threads; thread_id++) {
            if (thread_id > 0) {
                current_systems[thread_id]->set_prev_sys(current_systems[thread_id - 1]);
            } else {
                current_systems[thread_id]->set_prev_sys(last_system);
            }
        }

        // Parallel process system topology.
        std::vector<std::thread> threads;
        for (int thread_id = 0; thread_id < num_threads; thread_id++) {
            if (thread_id >= current_systems.size() || current_systems[thread_id]->atoms.size() == 0) {
                continue;
            }
            threads.emplace_back(&System::process_this, current_systems[thread_id]);
        }

        for (auto &thread : threads) {
            if (thread.joinable()) {
                thread.join();
            }
        }

        // Parallel process reaction issues

        threads.clear();
        for (int thread_id = 0; thread_id < num_threads; thread_id++) {
            if (thread_id >= current_systems.size() || current_systems[thread_id]->atoms.size() == 0) {
                continue;
            }

            threads.emplace_back(&System::process_reax, current_systems[thread_id]);
        }

        for (auto &thread : threads) {
            if (thread.joinable()) {
                thread.join();
            }
        }

        // Other issues, serial
        for (int thread_id = 0; thread_id < num_threads; thread_id++) {
            if (thread_id >= current_systems.size() || current_systems[thread_id]->atoms.size() == 0) {
                continue;
            }

            // current_systems[thread_id]->process_reax();

            if (if_dump_lammps_data) {
                std::string lammps_data_file =
                    output_dir + "frame_" + std::to_string(current_systems[thread_id]->frame_id) + ".data";
                current_systems[thread_id]->dump_lammps_data(lammps_data_file);
            }

            if (current_systems[thread_id]->frame_id == 1) {
                fmt::print("Atom Types: ");
                for (auto &pair : current_systems[thread_id]->type_itos) {
                    fmt::print("{}: {}, ", pair.first, pair.second);
                }
                fmt::print("\n");
                fmt::print("\n");

                is_first_frame = true;
            } else {
                is_first_frame = false;
            }

            current_systems[thread_id]->dump_bond_count(bond_count_filepath, is_first_frame);
            current_systems[thread_id]->dump_ring_count(ring_count_filepath, is_first_frame);

            current_systems[thread_id]->finish();
        }
    }
    fmt::print("\n\n");
    reax_species->analyze_frame_formulas();
}
