#pragma once

#include <algorithm>
#include <climits>
#include <cstdint>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "system.h"

// Forward declarations
struct Molecule;
class System;

struct TrackedMolecule {
    static int next_id;

    int id = 0;
    int start_frame = 0;
    int end_frame = INT_MAX;
    std::unordered_set<int> atom_ids;
    std::string formula;
    std::string smiles;          // empty for formula-level fallback
    std::uint64_t structure_id = 0;

    std::vector<TrackedMolecule*> predecessors;
    std::vector<TrackedMolecule*> successors;

    bool is_stable = false;
    bool is_processed = false;

    TrackedMolecule(int frame, const Molecule* mol);

    bool is_active() const { return end_frame == INT_MAX; }
    int lifetime(int current_frame) const {
        return is_active() ? current_frame - start_frame : end_frame - start_frame;
    }
    std::string to_string() const;
};

struct ReactionEvent {
    int event_id = 0;
    int frame_id = 0;
    float time_fs = 0.0f;

    std::vector<TrackedMolecule*> reactants;
    std::vector<TrackedMolecule*> products;

    ReactionEvent(int eid, int fid, float time);

    std::string get_reactants_string() const;  // formula list
    std::string get_products_string() const;
    bool check_atom_conservation() const;
};

class ReactionTracker {
private:
    std::unordered_map<int, TrackedMolecule*> atom_to_molecule;
    std::vector<TrackedMolecule*> all_molecules;
    std::unordered_set<TrackedMolecule*> active_molecules;
    std::vector<ReactionEvent> events;

    int stable_time_frames = 3;
    float timestep_fs = 0.25f;
    int sampling_frequency = 10;

    int next_event_id = 0;
    int last_processed_frame = 0;
    int transfer_edge_count = 0;

    TrackedMolecule* find_or_create_tracked(const Molecule* mol, int frame_id);
    void create_reaction_event(int frame_id,
                               const std::vector<TrackedMolecule*>& reactants,
                               const std::vector<TrackedMolecule*>& products);
    void detect_reactions_from_changes(
        int frame_id,
        const std::unordered_map<int, Molecule*>& current_atom_to_mol,
        const std::unordered_set<int>& changed_atoms);

public:
    ReactionTracker(int stable_frames = 3, float ts = 0.25f, int freq = 10);
    ~ReactionTracker();

    void process_frame(int frame_id, const std::vector<Molecule*>& molecules);
    void finalize(int last_frame);

    void save_raw_events(const std::string& filepath);
    void save_raw_event_pairs(const std::string& filepath);
    void save_transfer_flow(const std::string& filepath);

    int get_event_count() const { return static_cast<int>(events.size()); }
    int get_molecule_count() const { return static_cast<int>(all_molecules.size()); }
    int get_transfer_edge_count() const { return transfer_edge_count; }
    void brief_report() const;
};

extern bool FLAG_TRACK_REACTIONS;
extern int STABLE_TIME_FRAMES;
extern float TIMESTEP_FS;
extern int SAMPLING_FREQ;
