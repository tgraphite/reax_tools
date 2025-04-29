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
#include "molecule.h"
#include "reax_flow.h"
#include "reax_species.h"
#include "system.h"

class Universe {
   public:
    System *last_system = nullptr;
    std::vector<System *> current_systems;

    ReaxSpecies *reax_species = new ReaxSpecies();
    ReaxFlow *reax_flow = new ReaxFlow();

    Universe();
    ~Universe();

    void flush();

    void process_traj(std::string &file_path, std::vector<std::string> &type_names, const float &rvdw_scale,
                      const int &num_threads, const bool &if_dump_lammps_data);
};
