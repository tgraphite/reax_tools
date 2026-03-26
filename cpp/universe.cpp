#ifndef WASM_MODE
#include <pthread.h>
#endif

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "argparser.h"
#include "fmt/core.h"
#include "reax_counter.h"
#include "string_tools.h"
#include "universe.h"

Universe::Universe() {
    // FLAG_TRACK_REACTIONS is true by default
    reaction_tracker = new ReactionTracker(STABLE_TIME_FRAMES, TIMESTEP_FS, SAMPLING_FREQ);
}

#ifndef WASM_MODE
Universe::~Universe() {
    if (species_counter != nullptr) {
        delete species_counter;
        species_counter = nullptr;
    }

    if (reax_flow != nullptr) {
        delete reax_flow;
        reax_flow = nullptr;
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

    if (hash_counter != nullptr) {
        delete hash_counter;
        hash_counter = nullptr;
    }

    if (reaction_tracker != nullptr) {
        delete reaction_tracker;
        reaction_tracker = nullptr;
    }
}

void Universe::flush() {}
#else
Universe::~Universe() {
    if (system != nullptr) {
        delete system;
        system = nullptr;
    }

    if (last_system != nullptr) {
        delete last_system;
        last_system = nullptr;
    }

    if (species_counter != nullptr) {
        delete species_counter;
        species_counter = nullptr;
    }

    if (reax_flow != nullptr) {
        delete reax_flow;
        reax_flow = nullptr;
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

    if (hash_counter != nullptr) {
        delete hash_counter;
        hash_counter = nullptr;
    }

    if (reaction_tracker != nullptr) {
        delete reaction_tracker;
        reaction_tracker = nullptr;
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
#endif

#ifndef WSAM_MODE
template <typename T, typename Func>
void parallel_for_each(std::vector<T*>& objects, Func func) {
    size_t n = objects.size();
    std::vector<pthread_t> threads(n);

    struct ThreadCallArg {
        T* obj;
        Func func;
    };

    std::vector<ThreadCallArg> args(n);

    auto thread_func = [](void* arg) -> void* {
        ThreadCallArg* a = static_cast<ThreadCallArg*>(arg);
        if (a->obj) {
            (a->obj->*(a->func))();
        }
        return nullptr;
        };

    for (size_t i = 0; i < n; ++i) {
        args[i] = { objects[i], func };
        if (objects[i]) pthread_create(&threads[i], nullptr, thread_func, &args[i]);
    }
    for (size_t i = 0; i < n; ++i) {
        if (objects[i]) pthread_join(threads[i], nullptr);
    }
}
#endif

#ifndef WASM_MODE
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
            curr_system->reax_flow = reax_flow;
            curr_system->set_counters(this->species_counter, this->bond_counter, this->ring_counter,
                this->atom_bonded_num_counter, this->hash_counter);

            if (ends_with(INPUT_FILE, ".lammpstrj"))
                curr_system->load_lammpstrj(input_file);
            else if (ends_with(INPUT_FILE, ".xyz"))
                curr_system->load_xyz(input_file);

            if (curr_system->atoms.size() == 0) {
                delete curr_system;
                continue;
            };

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
        if (!FLAG_NO_REACTIONS) {
            parallel_for_each<System>(systems_to_process, &System::process_reax_flow);
        }

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

                std::unordered_set<unsigned int> initial_mol_hashes;
                for (auto& mol : curr_system->molecules) {
                    initial_mol_hashes.insert(mol->hash);
                }
                reax_flow->import_molecules(true, initial_mol_hashes);
            }

            if (curr_system->is_last_frame) {
                std::unordered_set<unsigned int> final_mol_hashes;
                for (auto& mol : curr_system->molecules) {
                    final_mol_hashes.insert(mol->hash);
                }
                reax_flow->import_molecules(false, final_mol_hashes);
            }

            curr_system->finish();
        }

        systems_to_process.clear();
    }
    fmt::print("\n\n");

    species_counter->analyze_frame_formulas();
    species_counter->brief_report();
    species_counter->save_file();

    reax_flow->save_graph();

    // after all batch, free all systems.
    for (auto it = frameid_system.begin(); it != frameid_system.end();) {
        delete it->second;              // Delete the instance
        it = frameid_system.erase(it);  // Remove the element from the map and advance iterator
    }

    // Finalize reaction tracking
    if (FLAG_TRACK_REACTIONS && reaction_tracker != nullptr) {
        reaction_tracker->finalize(curr_frame_id - 1);
        reaction_tracker->brief_report();
        reaction_tracker->save_events("reaction_track_events.csv");
        reaction_tracker->save_molecule_lifetimes("reaction_track_molecules.csv");
    }
}
#else
void Universe::process_traj() {
    int curr_frame_id = 1;
    int max_neigh = 10;

    std::ifstream file(INPUT_FILE);
    // std::string bond_count_filepath = output_dir + "bond_count.csv";
    // std::string ring_count_filepath = output_dir + "ring_count.csv";
    // std::string atom_bonded_num_count_filepath = output_dir + "atom_bonded_num_count.csv";

    // The highest calling stack, only do this once.
    while (file.is_open() && !file.eof()) {
        if (curr_frame_id > 1) {
            flush();
        }

        system = new System();
        system->frame_id = curr_frame_id;
        system->reax_flow = reax_flow;
        system->set_counters(this->species_counter, this->bond_counter, this->ring_counter,
            this->atom_bonded_num_counter, this->hash_counter);

        if (curr_frame_id == 1) {
            system->is_first_frame = true;
        }
        if (file.eof()) {
            system->is_last_frame = true;
        }

        if (ends_with(INPUT_FILE, ".lammpstrj"))
            system->load_lammpstrj(file);
        else if (ends_with(INPUT_FILE, ".xyz"))
            system->load_xyz(file);
        else {
            std::cerr << "Using unsupported format, check your file and suffix!";
            exit(1);
        }

        if (system->atoms.size() == 0) {
            continue;
        };
        curr_frame_id++;

        system->prev_sys = last_system;
        system->process_this();
        system->process_counters();
        system->process_reax_flow();

        if (FLAG_DUMP_STRUCTURE) {
            system->dump_lammps_data();
        }

        if (system->frame_id == 1) {
            fmt::print("Atom Types: ");
            for (auto& pair : system->type_itos) {
                fmt::print("{}: {}, ", pair.first, pair.second);
            }
            fmt::print("\n");
            fmt::print("\n");
        }

        system->finish();
    }
    fmt::print("\n\n");
    species_counter->analyze_frame_formulas();
}
#endif
