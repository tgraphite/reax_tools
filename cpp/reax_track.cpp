#include "reax_track.h"

#include <climits>
#include <queue>

#include "fmt/format.h"
#include "string_tools.h"

// External global variables (defined in argparser.cpp)
extern bool FLAG_TRACK_REACTIONS;
extern int STABLE_TIME_FRAMES;
extern float TIMESTEP_FS;
extern int SAMPLING_FREQ;

// Static member initialization
int TrackedMolecule::next_id = 0;

//=============================================================================
// TrackedMolecule Implementation
//=============================================================================

TrackedMolecule::TrackedMolecule(int frame, const Molecule* mol)
    : id(next_id++),
      start_frame(frame),
      end_frame(INT_MAX),
      formula(mol->formula),
      hash(mol->hash),
      is_stable(false),
      is_processed(false) {
    
    // Copy atom IDs
    for (int atom_id : mol->atom_ids) {
        atom_ids.insert(atom_id);
    }
}

std::string TrackedMolecule::to_string() const {
    return fmt::format("{}(id={},frames={}-{})", formula, id, start_frame, 
                       is_active() ? "inf" : std::to_string(end_frame));
}

//=============================================================================
// ReactionEvent Implementation
//=============================================================================

ReactionEvent::ReactionEvent(int eid, int fid, float time)
    : event_id(eid),
      frame_id(fid),
      time_fs(time),
      atom_transfer(0) {}

std::string ReactionEvent::get_reactants_string() const {
    std::string result;
    for (size_t i = 0; i < reactants.size(); ++i) {
        if (i > 0) result += "+";
        result += reactants[i]->formula;
    }
    return result;
}

std::string ReactionEvent::get_products_string() const {
    std::string result;
    for (size_t i = 0; i < products.size(); ++i) {
        if (i > 0) result += "+";
        result += products[i]->formula;
    }
    return result;
}

bool ReactionEvent::check_atom_conservation() const {
    std::unordered_set<int> reactant_atoms;
    std::unordered_set<int> product_atoms;
    
    for (const auto* mol : reactants) {
        reactant_atoms.insert(mol->atom_ids.begin(), mol->atom_ids.end());
    }
    for (const auto* mol : products) {
        product_atoms.insert(mol->atom_ids.begin(), mol->atom_ids.end());
    }
    
    return reactant_atoms == product_atoms;
}

//=============================================================================
// ReactionTracker Implementation
//=============================================================================

ReactionTracker::ReactionTracker(int stable_frames, float ts, int freq)
    : stable_time_frames(stable_frames),
      timestep_fs(ts),
      sampling_frequency(freq),
      next_event_id(0),
      last_processed_frame(0) {}

ReactionTracker::~ReactionTracker() {
    // Clean up all tracked molecules
    for (auto* mol : all_molecules) {
        delete mol;
    }
    all_molecules.clear();
}

TrackedMolecule* ReactionTracker::find_or_create_tracked(const Molecule* mol, int frame_id) {
    // First, check if any active molecule has the same set of atoms
    for (auto* tracked : active_molecules) {
        if (tracked->atom_ids == mol->atom_ids) {
            return tracked;
        }
    }
    
    // Create new tracked molecule
    TrackedMolecule* tracked = new TrackedMolecule(frame_id, mol);
    all_molecules.push_back(tracked);
    active_molecules.insert(tracked);
    
    return tracked;
}

bool ReactionTracker::should_mark_stable(TrackedMolecule* mol, int current_frame) {
    if (mol->is_stable) return false;
    if (mol->is_active()) {
        return (current_frame - mol->start_frame) >= stable_time_frames;
    }
    return (mol->end_frame - mol->start_frame) >= stable_time_frames;
}

void ReactionTracker::create_reaction_event(int frame_id,
                                           const std::vector<TrackedMolecule*>& reactants,
                                           const std::vector<TrackedMolecule*>& products) {
    // Calculate time in fs
    float time = frame_id * timestep_fs * sampling_frequency;
    
    ReactionEvent event(next_event_id++, frame_id, time);
    event.reactants = reactants;
    event.products = products;
    
    // Calculate atom transfer (intersection of atoms)
    std::unordered_set<int> reactant_atoms;
    std::unordered_set<int> product_atoms;
    
    for (const auto* mol : reactants) {
        reactant_atoms.insert(mol->atom_ids.begin(), mol->atom_ids.end());
    }
    for (const auto* mol : products) {
        product_atoms.insert(mol->atom_ids.begin(), mol->atom_ids.end());
    }
    
    // Atom transfer = size of intersection
    int common = 0;
    for (int atom : reactant_atoms) {
        if (product_atoms.count(atom)) common++;
    }
    event.atom_transfer = common;
    
    // Verify atom conservation
    if (!event.check_atom_conservation()) {
        fmt::print("Warning: Reaction event {} does not conserve atoms!\n", event.event_id);
    }
    
    events.push_back(std::move(event));
}

void ReactionTracker::detect_reactions_from_changes(
        int frame_id,
        const std::unordered_map<int, Molecule*>& current_atom_to_mol,
        const std::unordered_set<int>& changed_atoms) {
    
    std::unordered_set<int> processed_atoms;
    
    for (int start_atom : changed_atoms) {
        if (processed_atoms.count(start_atom)) continue;
        
        // Symmetric difference algorithm to find complete reaction cluster
        std::unordered_set<int> sym_diff = {start_atom};
        std::vector<TrackedMolecule*> reactants;
        std::vector<TrackedMolecule*> products;
        
        std::unordered_set<TrackedMolecule*> used_reactants;
        std::unordered_set<TrackedMolecule*> used_products;
        
        while (!sym_diff.empty()) {
            int atom_id = *sym_diff.begin();
            sym_diff.erase(sym_diff.begin());
            
            if (processed_atoms.count(atom_id)) continue;
            processed_atoms.insert(atom_id);
            
            // Check old molecule (from previous frame)
            auto old_it = atom_to_molecule.find(atom_id);
            if (old_it != atom_to_molecule.end()) {
                TrackedMolecule* old_tracked = old_it->second;
                if (!used_reactants.count(old_tracked)) {
                    reactants.push_back(old_tracked);
                    used_reactants.insert(old_tracked);
                    
                    // Add all atoms from this molecule to sym_diff
                    for (int a : old_tracked->atom_ids) {
                        if (!processed_atoms.count(a)) {
                            sym_diff.insert(a);
                        }
                    }
                }
            }
            
            // Check new molecule (in current frame)
            auto new_it = current_atom_to_mol.find(atom_id);
            if (new_it != current_atom_to_mol.end()) {
                Molecule* new_mol = new_it->second;
                TrackedMolecule* new_tracked = find_or_create_tracked(new_mol, frame_id);
                
                if (!used_products.count(new_tracked)) {
                    products.push_back(new_tracked);
                    used_products.insert(new_tracked);
                    
                    // Add all atoms from this molecule to sym_diff
                    for (int a : new_tracked->atom_ids) {
                        if (!processed_atoms.count(a)) {
                            sym_diff.insert(a);
                        }
                    }
                }
            }
        }
        
        // Only record if we have both reactants and products
        if (!reactants.empty() && !products.empty()) {
            // Set up predecessor/successor relationships
            for (auto* r : reactants) {
                for (auto* p : products) {
                    r->successors.push_back(p);
                    p->predecessors.push_back(r);
                }
            }
            
            // Mark reactants as ended
            for (auto* r : reactants) {
                if (r->is_active()) {
                    r->end_frame = frame_id;
                }
            }
            
            // Create reaction event
            create_reaction_event(frame_id, reactants, products);
        }
    }
}

void ReactionTracker::process_frame(int frame_id, const std::vector<Molecule*>& molecules) {
    if (!FLAG_TRACK_REACTIONS) return;
    
    // Build atom -> molecule mapping for current frame
    std::unordered_map<int, Molecule*> current_atom_to_mol;
    for (Molecule* mol : molecules) {
        for (int atom_id : mol->atom_ids) {
            current_atom_to_mol[atom_id] = mol;
        }
    }
    
    // First frame: just initialize mappings
    if (frame_id == 1 || atom_to_molecule.empty()) {
        for (Molecule* mol : molecules) {
            TrackedMolecule* tracked = find_or_create_tracked(mol, frame_id);
            for (int atom_id : mol->atom_ids) {
                atom_to_molecule[atom_id] = tracked;
            }
        }
        last_processed_frame = frame_id;
        return;
    }
    
    // Find atoms that changed molecular affiliation
    std::unordered_set<int> changed_atoms;
    
    for (const auto& [atom_id, new_mol] : current_atom_to_mol) {
        auto old_it = atom_to_molecule.find(atom_id);
        if (old_it != atom_to_molecule.end()) {
            TrackedMolecule* old_tracked = old_it->second;
            // Check if molecule hash changed
            if (old_tracked->hash != new_mol->hash) {
                changed_atoms.insert(atom_id);
            }
        } else {
            // New atom appeared (shouldn't happen in standard MD)
            changed_atoms.insert(atom_id);
        }
    }
    
    // Detect reactions using symmetric difference algorithm
    if (!changed_atoms.empty()) {
        detect_reactions_from_changes(frame_id, current_atom_to_mol, changed_atoms);
    }
    
    // Update atom_to_molecule for next frame
    atom_to_molecule.clear();
    for (Molecule* mol : molecules) {
        TrackedMolecule* tracked = find_or_create_tracked(mol, frame_id);
        for (int atom_id : mol->atom_ids) {
            atom_to_molecule[atom_id] = tracked;
        }
    }
    
    // Mark molecules that are no longer present as ended
    std::unordered_set<TrackedMolecule*> current_mols;
    for (Molecule* mol : molecules) {
        current_mols.insert(find_or_create_tracked(mol, frame_id));
    }
    
    for (auto* tracked : active_molecules) {
        if (!current_mols.count(tracked) && tracked->is_active()) {
            tracked->end_frame = frame_id;
        }
    }
    
    // Update active molecules set
    active_molecules = std::move(current_mols);
    last_processed_frame = frame_id;
}

void ReactionTracker::finalize(int last_frame) {
    if (!FLAG_TRACK_REACTIONS) return;
    
    // Mark all remaining active molecules as ended
    for (auto* mol : active_molecules) {
        if (mol->is_active()) {
            mol->end_frame = last_frame;
        }
    }
}

void ReactionTracker::save_events(const std::string& filepath) {
    if (events.empty()) return;
    
    // Count reaction frequencies and collect unique reactions
    struct ReactionRecord {
        std::string reactants_str;
        std::string products_str;
        int count;
        int first_frame;
        int last_frame;
    };
    
    std::map<std::string, ReactionRecord> reaction_map;
    
    for (const auto& event : events) {
        std::string key = event.get_reactants_string() + " -> " + event.get_products_string();
        auto it = reaction_map.find(key);
        if (it == reaction_map.end()) {
            reaction_map[key] = {
                event.get_reactants_string(),
                event.get_products_string(),
                1,
                event.frame_id,
                event.frame_id
            };
        } else {
            it->second.count++;
            it->second.last_frame = event.frame_id;
        }
    }
    
    // Convert to vector and sort by frequency (descending)
    std::vector<ReactionRecord> sorted_reactions;
    for (const auto& [key, record] : reaction_map) {
        sorted_reactions.push_back(record);
    }
    std::sort(sorted_reactions.begin(), sorted_reactions.end(),
              [](const auto& a, const auto& b) { return a.count > b.count; });
    
    FILE* fp = create_file(filepath);
    
    // Header (without atom_transfer)
    fmt::print(fp, "rank,frequency,first_frame,last_frame,n_reactants,n_products,reactants,products\n");
    
    // Reactions sorted by frequency
    int rank = 1;
    for (const auto& rec : sorted_reactions) {
        // Count reactants and products
        int n_reactants = 1;
        int n_products = 1;
        for (char c : rec.reactants_str) if (c == '+') n_reactants++;
        for (char c : rec.products_str) if (c == '+') n_products++;
        
        fmt::print(fp, "{},{},{},{},{},{},\"{}\",\"{}\"\n",
                   rank++,
                   rec.count,
                   rec.first_frame,
                   rec.last_frame,
                   n_reactants,
                   n_products,
                   rec.reactants_str,
                   rec.products_str);
    }
    
    fclose(fp);
    fmt::print("Saved {} unique reaction types (sorted by frequency) to {}\n", sorted_reactions.size(), filepath);
}

void ReactionTracker::save_molecule_lifetimes(const std::string& filepath) {
    if (all_molecules.empty()) return;
    
    FILE* fp = create_file(filepath);
    
    // Header
    fmt::print(fp, "mol_id,formula,first_frame,last_frame,lifetime_frames,n_predecessors,n_successors\n");
    
    // Molecules
    for (const auto* mol : all_molecules) {
        int lifetime = mol->end_frame - mol->start_frame;
        fmt::print(fp, "{},{},{},{},{},{},{}\n",
                   mol->id,
                   mol->formula,
                   mol->start_frame,
                   mol->is_active() ? last_processed_frame : mol->end_frame,
                   lifetime,
                   mol->predecessors.size(),
                   mol->successors.size());
    }
    
    fclose(fp);
    fmt::print("Saved {} molecule lifetimes to {}\n", all_molecules.size(), filepath);
}

void ReactionTracker::brief_report() const {
    fmt::print("\n=== Reaction Tracker Report ===\n");
    fmt::print("Total molecules tracked: {}\n", all_molecules.size());
    fmt::print("Active molecules: {}\n", active_molecules.size());
    fmt::print("Reaction events detected: {}\n", events.size());
    
    if (!events.empty()) {
        fmt::print("\nTop 10 most common reactions:\n");
        
        // Count reaction types
        std::map<std::string, int> reaction_counts;
        for (const auto& event : events) {
            std::string key = event.get_reactants_string() + " -> " + event.get_products_string();
            reaction_counts[key]++;
        }
        
        // Sort by count
        std::vector<std::pair<std::string, int>> sorted(reaction_counts.begin(), reaction_counts.end());
        std::sort(sorted.begin(), sorted.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        
        int max_display = std::min(10, (int)sorted.size());
        for (int i = 0; i < max_display; ++i) {
            fmt::print("  {}: {} times\n", sorted[i].first, sorted[i].second);
        }
    }
}
