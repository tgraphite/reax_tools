#include "universe.h"

#include <omp.h>

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

void Universe::process_traj(std::string &file_path, std::vector<std::string> &type_names, const float &rvdw_scale,
                            const int &num_threads, const bool &if_dump_lammps_data) {
    std::ifstream file(file_path);
    int curr_frame = 1;
    float neigh_radius = 2.5 * rvdw_scale;

    while (file.is_open() and !file.eof()) {
        if (curr_frame > 1) {
            flush();
        }

        // Serial read. systems readed in the right order.
        for (int curr_thread = 0; curr_thread < num_threads; curr_thread++) {
            System *curr_sys = new System();
            curr_sys->set_types(type_names);
            current_systems.push_back(curr_sys);

            if (file_path.ends_with(".lammpstrj"))
                curr_sys->load_lammpstrj(file);
            else if (file_path.ends_with(".xyz"))
                curr_sys->load_xyz(file);

            if (curr_sys->atoms.size() == 0) continue;
        }

        // Parallel process system topology.
        std::vector<std::thread> threads;
        for (int thread_id = 0; thread_id < num_threads; thread_id++) {
            if (thread_id >= current_systems.size() || current_systems[thread_id]->atoms.size() == 0) {
                continue;
            }

            threads.push_back(std::thread([&, thread_id]() {
                System *curr_sys = current_systems[thread_id];
                curr_sys->search_neigh(neigh_radius, 10);
                curr_sys->build_bonds_by_radius(rvdw_scale);
                curr_sys->build_molecules();
            }));
        }

        //
        for (auto &thread : threads) {
            if (thread.joinable()) {
                thread.join();
            }
        }

        // Other issues
        for (int curr_thread = 0; curr_thread < num_threads; curr_thread++) {
            if (curr_thread >= current_systems.size() || current_systems[curr_thread]->atoms.size() == 0) {
                continue;
            }

            System *prev_sys;
            System *curr_sys = current_systems[curr_thread];

            std::vector<std::string> frame_formulas(curr_sys->molecules.size());
            for (size_t i = 0; i < curr_sys->molecules.size(); i++) {
                frame_formulas[i] = curr_sys->molecules[i]->formula;
            }
            reax_species->import_frame_formulas(frame_formulas);

            if (curr_thread == 0) {
                update_reax_flow(last_system, curr_sys, curr_frame);
            } else {
                update_reax_flow(current_systems[curr_thread - 1], curr_sys, curr_frame);
            }

            fmt::print("Frame: {} ", curr_frame);
            if (if_dump_lammps_data) {
                std::string lammps_data_file =
                    file_path.substr(0, file_path.find_last_of(".")) + "_" + std::to_string(curr_frame) + ".data";
                curr_sys->finish(lammps_data_file);
            } else {
                curr_sys->finish();
            }

            curr_frame++;
        }
    }
    reax_species->analyze_frame_formulas();
}

void Universe::update_reax_flow(System *prev_sys, System *curr_sys, const int &curr_frame) {
    if (prev_sys == nullptr || curr_sys == nullptr) {
        return;
    }

    for (auto &prev_mol : prev_sys->molecules) {
        if (prev_mol == nullptr) continue;

        // Ignore single atom molecule.
        if (prev_mol->atom_ids.size() == 1) continue;

        // Find the most similar molecule in current frame.
        Molecule *best_match = nullptr;
        float best_similarity = 0.0f;

        std::vector<int> intersection;
        std::vector<int> union_set;

        // Optimization: stop searching once a similarity > 0.5 is found.
        bool found_match = false;

        for (auto &curr_mol : curr_sys->molecules) {
            if (curr_mol == nullptr) continue;

            // Ignore single atom molecule.
            // if (curr_mol->atom_ids.size() == 1) continue;

            // If the formula is the same, consider it the same molecule, skip.
            if (prev_mol->formula == curr_mol->formula) continue;

            // Quick filter: skip if the size difference is too large.
            // if ((curr_mol->atom_ids.size() / prev_mol->atom_ids.size() >= 2) ||
            //     (prev_mol->atom_ids.size() / curr_mol->atom_ids.size() >= 2))
            //     continue;

            // Calculate similarity: intersection / union.
            intersection.clear();
            union_set.clear();

            std::set_intersection(prev_mol->atom_ids.begin(), prev_mol->atom_ids.end(), curr_mol->atom_ids.begin(),
                                  curr_mol->atom_ids.end(), back_inserter(intersection));

            // Quick check: if intersection is empty, similarity is 0.
            if (intersection.empty()) continue;

            // Quick check: if intersection size equals prev_mol size and
            // curr_mol size, then they are the same molecule.
            if (intersection.size() == prev_mol->atom_ids.size() && intersection.size() == curr_mol->atom_ids.size())
                continue;  // Same molecule, not reaction.

            std::set_union(prev_mol->atom_ids.begin(), prev_mol->atom_ids.end(), curr_mol->atom_ids.begin(),
                           curr_mol->atom_ids.end(), back_inserter(union_set));

            float similarity = float(intersection.size()) / float(union_set.size());

            // If similarity exceeds the threshold, record and stop searching.
            if (similarity > 0.5) {
                best_match = curr_mol;
                best_similarity = similarity;
                found_match = true;
                break;  // If a good match is found, stop searching.
            }

            // If similarity exceeds the threshold, record and stop searching.
            else if (similarity > best_similarity) {
                best_match = curr_mol;
                best_similarity = similarity;
            }
        }

        // If found a match and not the same molecule (similarity between 0.5
        // and 1.0), record the reaction.
        if (best_match && best_similarity >= 0.01 && best_similarity < 1.0) {
            reax_flow->add_reaction(curr_frame, prev_mol, best_match);
        }
    }
}
