#include <pthread.h>

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
    if (reax_species != nullptr) {
        delete reax_species;
        reax_species = nullptr;
    }

    if (reax_flow != nullptr) {
        delete reax_flow;
        reax_flow = nullptr;
    }
}

void Universe::flush() {}

template <typename T, typename Func>
void parallel_for_each(std::vector<T *> &objects, Func func) {
    size_t n = objects.size();
    std::vector<pthread_t> threads(n);

    struct ThreadCallArg {
        T *obj;
        Func func;
    };

    std::vector<ThreadCallArg> args(n);

    auto thread_func = [](void *arg) -> void * {
        ThreadCallArg *a = static_cast<ThreadCallArg *>(arg);
        if (a->obj) {
            (a->obj->*(a->func))();
        }
        return nullptr;
    };

    for (size_t i = 0; i < n; ++i) {
        args[i] = {objects[i], func};
        if (objects[i]) pthread_create(&threads[i], nullptr, thread_func, &args[i]);
    }
    for (size_t i = 0; i < n; ++i) {
        if (objects[i]) pthread_join(threads[i], nullptr);
    }
}

void Universe::process_traj(std::string &file_path, std::string &output_dir, std::vector<std::string> &type_names,
                            const float &rvdw_scale, const int &num_threads, const bool &if_dump_lammps_data) {
    float neigh_radius = 2.4 * rvdw_scale;
    int max_neigh = 10;
    int curr_frame_id = 1;
    bool is_first_frame = true;

    std::ifstream file(file_path);
    std::string bond_count_filepath = output_dir + "bond_count.csv";
    std::string ring_count_filepath = output_dir + "ring_count.csv";
    std::string atom_bonded_num_count_filepath = output_dir + "atom_bonded_num_count.csv";

    std::map<int, System *> frameid_system;
    std::vector<System *> systems_to_process;

    while (file.is_open() && !file.eof()) {
        // flush
        if (curr_frame_id > 1) {
            // Iterate safely and erase elements from the map if their instance should be destroyed.
            for (auto it = frameid_system.begin(); it != frameid_system.end();) {
                if (it->second->to_destroy) {
                    delete it->second;              // Delete the instance
                    it = frameid_system.erase(it);  // Remove the element from the map and advance iterator
                } else {
                    ++it;
                }
            }
        }

        // Initialize system and read data.
        for (size_t thread_id = 0; thread_id < num_threads; thread_id++) {
            System *curr_system = new System();
            curr_system->set_frame_id(curr_frame_id);
            curr_system->set_types(type_names);
            curr_system->set_max_neigh(max_neigh);
            curr_system->set_rvdw_scale(rvdw_scale);
            curr_system->set_reax_flow(this->reax_flow);
            curr_system->set_reax_species(this->reax_species);

            if (ends_with(file_path, ".lammpstrj"))
                curr_system->load_lammpstrj(file);
            else if (ends_with(file_path, ".xyz"))
                curr_system->load_xyz(file);

            if (curr_system->atoms.size() == 0) {
                delete curr_system;
                continue;
            };

            frameid_system[curr_frame_id] = curr_system;
            systems_to_process.push_back(curr_system);
            curr_frame_id++;
        }

        parallel_for_each<System>(systems_to_process, &System::process_this);
        for (auto &curr_system : systems_to_process) {
            int frame_id = curr_system->frame_id;
            if (frame_id == 1) curr_system->set_prev_sys(nullptr);
            // process_reax will skip compuations related to prev_sys.
            else
                curr_system->set_prev_sys(frameid_system[frame_id - 1]);
        }
        parallel_for_each<System>(systems_to_process, &System::process_reax);

        for (auto &curr_system : systems_to_process) {
            if (curr_system->prev_sys)  // prev_sys of the first frame is nullptr
                curr_system->prev_sys->to_destroy = true;

            if (if_dump_lammps_data) {
                std::string lammps_data_file = output_dir + "frame_" + std::to_string(curr_system->frame_id) + ".data";
                curr_system->dump_lammps_data(lammps_data_file);
            }

            if (curr_system->frame_id == 1) {
                fmt::print("Atom Types: ");
                for (auto &pair : curr_system->type_itos) {
                    fmt::print("{}: {}, ", pair.first, pair.second);
                }

                fmt::print("\n");

                fmt::print("Bond radius: ");
                for (auto &pair : curr_system->bond_radius_sq) {
                    fmt::print("{}-{} {:.3f}, ", curr_system->type_itos[pair.first.first],
                               curr_system->type_itos[pair.first.second], std::sqrt(pair.second));
                }

                fmt::print("\n");

                is_first_frame = true;
            } else {
                is_first_frame = false;
            }

            curr_system->dump_bond_count(bond_count_filepath, is_first_frame);
            curr_system->dump_ring_count(ring_count_filepath, is_first_frame);
            curr_system->dump_atom_bonded_num_count(atom_bonded_num_count_filepath, is_first_frame);
            curr_system->finish();
        }

        systems_to_process.clear();
    }
    fmt::print("\n\n");
    reax_species->analyze_frame_formulas();

    // after all batch, free all systems.
    for (auto it = frameid_system.begin(); it != frameid_system.end();) {
        delete it->second;              // Delete the instance
        it = frameid_system.erase(it);  // Remove the element from the map and advance iterator
    }
}
