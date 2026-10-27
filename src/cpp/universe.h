#pragma once
#include <cstddef>
#include <string>

#include "reax_counter.h"
#include "reax_track.h"
#include "system.h"

class Universe {
   public:
    System* last_system = nullptr;
    System* system = nullptr;

    SpeciesCounter* species_counter = new SpeciesCounter();
    ReactionTracker* reaction_tracker = nullptr;
    Counter<int>* bond_counter = new Counter<int>();
    Counter<int>* ring_counter = new Counter<int>();
    Counter<int>* atom_bonded_num_counter = new Counter<int>();

    int frames_processed = 0;
    std::size_t atoms_first_frame = 0;
    std::size_t atoms_last_frame = 0;

    Universe();
    ~Universe();

    void flush();
    void process_traj();
};
