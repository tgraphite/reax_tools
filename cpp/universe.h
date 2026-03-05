#pragma once
#include <string>

#include "reax_counter.h"
#include "reax_flow.h"
#include "system.h"

class Universe {
   public:
    System* last_system = nullptr;
    System* system = nullptr;

    SpeciesCounter* species_counter = new SpeciesCounter();
    ReaxFlow* reax_flow = new ReaxFlow();
    Counter<int>* bond_counter = new Counter<int>();
    Counter<int>* ring_counter = new Counter<int>();
    Counter<int>* atom_bonded_num_counter = new Counter<int>();
    Counter<int>* hash_counter = new Counter<int>();

    Universe();
    ~Universe();

    void flush();
    void process_traj();
};
