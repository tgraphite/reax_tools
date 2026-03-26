#pragma once

#include <algorithm>
#include <climits>
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

/**
 * @brief TrackedMolecule - A molecule with lifecycle tracking
 * 
 * Similar to ChemTrayzer's MDMolecule, this class tracks a molecule's
 * existence over multiple frames, including its predecessors and successors.
 */
struct TrackedMolecule {
    static int next_id;  // Global counter for unique IDs
    
    int id;                          // Unique molecule ID
    int start_frame;                 // First frame where this molecule appears
    int end_frame;                   // Last frame (INT_MAX if still active)
    std::unordered_set<int> atom_ids;  // Set of atom IDs in this molecule
    std::string formula;             // Chemical formula
    unsigned int hash;               // Hash of formula
    
    // Network connections (similar to ChemTrayzer's MDMoleculeNetwork)
    std::vector<TrackedMolecule*> predecessors;
    std::vector<TrackedMolecule*> successors;
    
    // Tracking state
    bool is_stable;                  // Whether molecule has reached stable lifetime
    bool is_processed;               // Whether reaction involving this mol has been recorded
    
    TrackedMolecule(int frame, const Molecule* mol);
    ~TrackedMolecule() = default;
    
    // Check if molecule is still active (not ended)
    bool is_active() const { return end_frame == INT_MAX; }
    
    // Get lifetime in frames
    int lifetime(int current_frame) const {
        if (is_active()) return current_frame - start_frame;
        return end_frame - start_frame;
    }
    
    // String representation for output
    std::string to_string() const;
};

/**
 * @brief ReactionEvent - A single reaction event detected between frames
 * 
 * Records the reactants, products, frame/time, and atom transfer information.
 */
struct ReactionEvent {
    int event_id;                    // Unique event ID
    int frame_id;                    // Frame where reaction was detected
    float time_fs;                   // Simulation time in femtoseconds
    
    std::vector<TrackedMolecule*> reactants;
    std::vector<TrackedMolecule*> products;
    
    int atom_transfer;               // Number of atoms transferred
    
    // For detailed analysis
    std::set<std::pair<int, int>> formed_bonds;
    std::set<std::pair<int, int>> broken_bonds;
    
    ReactionEvent(int eid, int fid, float time);
    
    // Get string representation of reactants/products
    std::string get_reactants_string() const;
    std::string get_products_string() const;
    
    // Verify atom conservation
    bool check_atom_conservation() const;
};

/**
 * @brief ReactionTracker - Main class for tracking molecular evolution and reactions
 * 
 * Implements the ChemTrayzer-like algorithm for detecting reactions from
 * molecular dynamics trajectories without bond order information.
 */
class ReactionTracker {
private:
    // Current state: atom ID -> TrackedMolecule
    std::unordered_map<int, TrackedMolecule*> atom_to_molecule;
    
    // All molecules ever tracked
    std::vector<TrackedMolecule*> all_molecules;
    
    // Currently active molecules (still existing in current frame)
    std::unordered_set<TrackedMolecule*> active_molecules;
    
    // Detected reaction events
    std::vector<ReactionEvent> events;
    
    // Configuration
    int stable_time_frames;          // TStable: minimum frames for molecule to be stable
    float timestep_fs;               // MD timestep in fs
    int sampling_frequency;          // Frame sampling frequency
    
    // Internal tracking
    int next_event_id;
    int last_processed_frame;
    
    // Helper: Find or create tracked molecule from System's Molecule
    TrackedMolecule* find_or_create_tracked(const Molecule* mol, int frame_id);
    
    // Helper: Check if a molecule should be marked as stable
    bool should_mark_stable(TrackedMolecule* mol, int current_frame);
    
    // Helper: Create reaction event from reactants and products
    void create_reaction_event(int frame_id, 
                               const std::vector<TrackedMolecule*>& reactants,
                               const std::vector<TrackedMolecule*>& products);
    
    // Helper: Detect reactions using symmetric difference algorithm
    void detect_reactions_from_changes(int frame_id, 
                                       const std::unordered_map<int, Molecule*>& current_atom_to_mol,
                                       const std::unordered_set<int>& changed_atoms);
    
    // Helper: Mark stable molecules and create events
    void process_stable_molecules(int current_frame);
    
public:
    ReactionTracker(int stable_frames = 3, float ts = 0.25f, int freq = 10);
    ~ReactionTracker();
    
    // Process a single frame - main entry point
    void process_frame(int frame_id, const std::vector<Molecule*>& molecules);
    
    // Finalize and save results
    void finalize(int last_frame);
    void save_events(const std::string& filepath);
    void save_molecule_lifetimes(const std::string& filepath);
    
    // Getters
    int get_event_count() const { return events.size(); }
    int get_molecule_count() const { return all_molecules.size(); }
    
    // Brief report
    void brief_report() const;
};

// Global pointer for easy access (similar to other globals in reax_tools)
extern bool FLAG_TRACK_REACTIONS;
extern int STABLE_TIME_FRAMES;
extern float TIMESTEP_FS;
extern int SAMPLING_FREQ;
