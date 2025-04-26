#pragma once
#include <algorithm>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#include "atom.h"
#include "kdtree.h"
#include "molecule.h"
#include "reax_flow.h"
#include "reax_species.h"
#include "system.h"

class Universe {
   public:
    std::vector<System *> current_systems;
    System *last_system = nullptr;

    ReaxSpecies *reax_species = new ReaxSpecies();
    ReaxFlow *reax_flow = new ReaxFlow();

    Universe();
    ~Universe();

    float compute_similarity(std::unordered_set<int> &prev_set, std::unordered_set<int> &curr_set);
    void update_reax_flow(System *prev_sys, System *curr_sys, const int &curr_frame);
    void flush();

    void process_traj(std::string &file_path, std::vector<std::string> &type_names, const float &rvdw_scale,
                      const int &num_threads, const bool &if_dump_lammps_data);

    void dump_final_frame_data(std::string &filepath);
};
