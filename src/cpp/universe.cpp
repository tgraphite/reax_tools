#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#ifdef REAX_ENABLE_THREADS
#include <thread>
#endif

#include "argparser.h"
#include "fmt/core.h"
#include "reax_counter.h"
#include "string_tools.h"
#include "universe.h"

Universe::Universe() {
    // FLAG_TRACK_REACTIONS is true by default
    reaction_tracker = new ReactionTracker(STABLE_TIME_FRAMES, TIMESTEP_FS, SAMPLING_FREQ);
}

Universe::~Universe() {
    if (species_counter != nullptr) {
        delete species_counter;
        species_counter = nullptr;
    }

    if (bond_counter != nullptr) {
        delete bond_counter;
        bond_counter = nullptr;
    }

    if (ring_counter != nullptr) {
        delete ring_counter;
        ring_counter = nullptr;
    }

    if (atom_bonded_num_counter != nullptr) {
        delete atom_bonded_num_counter;
        atom_bonded_num_counter = nullptr;
    }

    if (reaction_tracker != nullptr) {
        delete reaction_tracker;
        reaction_tracker = nullptr;
    }
}

void Universe::flush() {}

#ifdef REAX_ENABLE_THREADS
template <typename T, typename Func>
void parallel_for_each(std::vector<T*>& objects, Func func) {
    std::vector<std::thread> threads;
    threads.reserve(objects.size());
    for (T* obj : objects) {
        if (obj) threads.emplace_back([obj, func]() { (obj->*func)(); });
    }
    for (auto& thread : threads) {
        if (thread.joinable()) thread.join();
    }
}
#else
template <typename T, typename Func>
void parallel_for_each(std::vector<T*>& objects, Func func) {
    for (T* obj : objects) {
        if (obj) (obj->*func)();
    }
}
#endif

void Universe::process_traj() {
    int max_neigh = 10;
    int curr_frame_id = 1;

    std::ifstream input_file(INPUT_FILE);

    std::map<int, System*> frameid_system;
    std::vector<System*> systems_to_process;

    while (input_file.is_open() && !input_file.eof()) {
        // flush
        if (curr_frame_id > 1) {
            // Iterate safely and erase elements from the map if their instance should be destroyed.
            for (auto it = frameid_system.begin(); it != frameid_system.end();) {
                if (it->second->to_destroy) {
                    delete it->second;              // Delete the instance
                    it = frameid_system.erase(it);  // Remove the element from the map and advance iterator
                }
                else {
                    ++it;
                }
            }
        }

        // Initialize system and read data.
        for (size_t thread_id = 0; thread_id < NUM_THREADS; thread_id++) {
            System* curr_system = new System();

            curr_system->frame_id = curr_frame_id;
            curr_system->set_counters(this->species_counter, this->bond_counter, this->ring_counter,
                this->atom_bonded_num_counter);

            if (ends_with(INPUT_FILE, ".lammpstrj"))
                curr_system->load_lammpstrj(input_file);
            else if (ends_with(INPUT_FILE, ".xyz"))
                curr_system->load_xyz(input_file);

            if (curr_system->atoms.size() == 0) {
                delete curr_system;
                continue;
            };

            if (curr_frame_id == 1) {
                atoms_first_frame = curr_system->atoms.size();
            }
            atoms_last_frame = curr_system->atoms.size();

            if (curr_frame_id == 1) {
                curr_system->is_first_frame = true;
            }
            else {
                curr_system->is_last_frame = true;
            }

            frameid_system[curr_frame_id] = curr_system;
            systems_to_process.push_back(curr_system);
            curr_frame_id++;
        }

        parallel_for_each<System>(systems_to_process, &System::process_this);
        for (auto& curr_system : systems_to_process) {
            int frame_id = curr_system->frame_id;
            if (frame_id == 1) curr_system->prev_sys = nullptr;
            // process_reax will skip compuations related to prev_sys.
            else {
                curr_system->prev_sys = frameid_system[frame_id - 1];
                curr_system->prev_sys->is_last_frame = false;  // if a system can be prev_sys, it is not the last frame.
            }
        }
        parallel_for_each<System>(systems_to_process, &System::process_counters);

        // Process reaction tracking (must be sequential by frame)
        if (FLAG_TRACK_REACTIONS && reaction_tracker != nullptr) {
            // Sort systems by frame_id to ensure correct order
            std::sort(systems_to_process.begin(), systems_to_process.end(),
                [](System* a, System* b) { return a->frame_id < b->frame_id; });
            for (auto& curr_system : systems_to_process) {
                reaction_tracker->process_frame(curr_system->frame_id, curr_system->molecules);
            }
        }

        for (auto& curr_system : systems_to_process) {
            if (curr_system->prev_sys)  // prev_sys of the first frame is nullptr
                curr_system->prev_sys->to_destroy = true;

            if (FLAG_DUMP_STRUCTURE) {
                curr_system->dump_lammps_data();
            }

            if (curr_system->is_first_frame) {
                fmt::print("Atom Types: ");
                for (auto& pair : curr_system->type_itos) {
                    fmt::print("{}: {}, ", pair.first, pair.second);
                }

                fmt::print("\n");

                fmt::print("Bond radius: ");
                for (auto& pair : curr_system->bond_radius_sq) {
                    fmt::print("{}-{} {:.3f}, ", curr_system->type_itos[pair.first.first],
                        curr_system->type_itos[pair.first.second], std::sqrt(pair.second));
                }

                fmt::print("\n");
            }

            curr_system->finish();
        }

        systems_to_process.clear();
    }
    frames_processed = curr_frame_id - 1;
    fmt::print("\n\n");

    species_counter->analyze_frame_species();
    species_counter->save_species_count();
    species_counter->save_molecules();

    // after all batch, free all systems.
    for (auto it = frameid_system.begin(); it != frameid_system.end();) {
        delete it->second;              // Delete the instance
        it = frameid_system.erase(it);  // Remove the element from the map and advance iterator
    }

    // Finalize reaction tracking
    if (FLAG_TRACK_REACTIONS && reaction_tracker != nullptr) {
        reaction_tracker->finalize(curr_frame_id - 1);
        reaction_tracker->save_raw_events("reaction_events.csv");
        reaction_tracker->save_raw_event_pairs("reaction_event_pairs.csv");
        reaction_tracker->save_transfer_flow("transfer_flow.csv");
    }
}
