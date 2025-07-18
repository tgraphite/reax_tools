#pragma once
#include <string>

#include "reax_counter.h"
#include "reax_flow.h"
#include "system.h"

class Universe {
  public:
    System* last_system = nullptr;
    System* system = nullptr;

    ReaxCounter* reax_counter = new ReaxCounter();
    ReaxFlow* reax_flow = new ReaxFlow();
    Counter<int>* bond_counter = new Counter<int>();
    Counter<int>* ring_counter = new Counter<int>();
    Counter<int>* atom_bonded_num_counter = new Counter<int>();

    Universe();
    ~Universe();

    void flush();

    void process_traj(std::string& file_path, std::string& output_dir, std::vector<std::string>& type_names,
                      const float& rvdw_scale, const int& num_threads, const bool& if_dump_lammps_data,
                      const int& dump_data_frame_step, const bool& if_mark_ring_atoms, const bool& if_no_reax_flow);
};
