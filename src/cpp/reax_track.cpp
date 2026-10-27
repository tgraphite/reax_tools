#include "reax_track.h"

#include <climits>
#include <map>
#include <string>
#include <vector>

#include "argparser.h"
#include "fmt/format.h"
#include "string_tools.h"

extern bool FLAG_TRACK_REACTIONS;
extern int STABLE_TIME_FRAMES;
extern float TIMESTEP_FS;
extern int SAMPLING_FREQ;

namespace {

std::string join_ids(const std::vector<TrackedMolecule*>& molecules) {
    std::string out;
    for (std::size_t i = 0; i < molecules.size(); ++i) {
        if (i > 0) out += ";";
        out += std::to_string(molecules[i]->structure_id);
    }
    return out;
}

std::string join_formulas(const std::vector<TrackedMolecule*>& molecules) {
    std::string out;
    for (std::size_t i = 0; i < molecules.size(); ++i) {
        if (i > 0) out += ";";
        out += molecules[i]->formula;
    }
    return out;
}

std::string join_smiles(const std::vector<TrackedMolecule*>& molecules) {
    std::string out;
    for (std::size_t i = 0; i < molecules.size(); ++i) {
        if (i > 0) out += ";";
        out += molecules[i]->smiles;
    }
    return out;
}

}  // namespace

int TrackedMolecule::next_id = 0;

TrackedMolecule::TrackedMolecule(int frame, const Molecule* mol)
    : id(next_id++),
      start_frame(frame),
      end_frame(INT_MAX),
      formula(mol->formula),
      smiles(mol->smiles),
      structure_id(mol->structure_id) {
    atom_ids = mol->atom_ids;
}

std::string TrackedMolecule::to_string() const {
    return fmt::format("{}(id={},frames={}-{})", formula, id, start_frame,
                       is_active() ? "inf" : std::to_string(end_frame));
}

ReactionEvent::ReactionEvent(int eid, int fid, float time)
    : event_id(eid), frame_id(fid), time_fs(time) {}

std::string ReactionEvent::get_reactants_string() const {
    std::string result;
    for (std::size_t i = 0; i < reactants.size(); ++i) {
        if (i > 0) result += "+";
        result += reactants[i]->formula;
    }
    return result;
}

std::string ReactionEvent::get_products_string() const {
    std::string result;
    for (std::size_t i = 0; i < products.size(); ++i) {
        if (i > 0) result += "+";
        result += products[i]->formula;
    }
    return result;
}

bool ReactionEvent::check_atom_conservation() const {
    std::unordered_set<int> reactant_atoms;
    std::unordered_set<int> product_atoms;
    for (const auto* mol : reactants) reactant_atoms.insert(mol->atom_ids.begin(), mol->atom_ids.end());
    for (const auto* mol : products) product_atoms.insert(mol->atom_ids.begin(), mol->atom_ids.end());
    return reactant_atoms == product_atoms;
}

ReactionTracker::ReactionTracker(int stable_frames, float ts, int freq)
    : stable_time_frames(stable_frames),
      timestep_fs(ts),
      sampling_frequency(freq),
      next_event_id(0),
      last_processed_frame(0) {}

ReactionTracker::~ReactionTracker() {
    for (auto* mol : all_molecules) delete mol;
    all_molecules.clear();
}

TrackedMolecule* ReactionTracker::find_or_create_tracked(const Molecule* mol, int frame_id) {
    for (auto* tracked : active_molecules) {
        if (tracked->atom_ids == mol->atom_ids && tracked->structure_id == mol->structure_id) {
            return tracked;
        }
    }

    TrackedMolecule* tracked = new TrackedMolecule(frame_id, mol);
    all_molecules.push_back(tracked);
    active_molecules.insert(tracked);
    return tracked;
}

void ReactionTracker::create_reaction_event(
    int frame_id,
    const std::vector<TrackedMolecule*>& reactants,
    const std::vector<TrackedMolecule*>& products) {
    float time = frame_id * timestep_fs * sampling_frequency;
    ReactionEvent event(next_event_id++, frame_id, time);
    event.reactants = reactants;
    event.products = products;
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

            auto old_it = atom_to_molecule.find(atom_id);
            if (old_it != atom_to_molecule.end()) {
                TrackedMolecule* old_tracked = old_it->second;
                if (!used_reactants.count(old_tracked)) {
                    reactants.push_back(old_tracked);
                    used_reactants.insert(old_tracked);
                    for (int a : old_tracked->atom_ids) {
                        if (!processed_atoms.count(a)) sym_diff.insert(a);
                    }
                }
            }

            auto new_it = current_atom_to_mol.find(atom_id);
            if (new_it != current_atom_to_mol.end()) {
                TrackedMolecule* new_tracked = find_or_create_tracked(new_it->second, frame_id);
                if (!used_products.count(new_tracked)) {
                    products.push_back(new_tracked);
                    used_products.insert(new_tracked);
                    for (int a : new_tracked->atom_ids) {
                        if (!processed_atoms.count(a)) sym_diff.insert(a);
                    }
                }
            }
        }

        if (!reactants.empty() && !products.empty()) {
            for (auto* r : reactants) {
                for (auto* p : products) {
                    r->successors.push_back(p);
                    p->predecessors.push_back(r);
                }
            }
            for (auto* r : reactants) {
                if (r->is_active()) r->end_frame = frame_id;
            }
            create_reaction_event(frame_id, reactants, products);
        }
    }
}

void ReactionTracker::process_frame(int frame_id, const std::vector<Molecule*>& molecules) {
    if (!FLAG_TRACK_REACTIONS) return;

    std::unordered_map<int, Molecule*> current_atom_to_mol;
    for (Molecule* mol : molecules) {
        for (int atom_id : mol->atom_ids) current_atom_to_mol[atom_id] = mol;
    }

    if (frame_id == 1 || atom_to_molecule.empty()) {
        for (Molecule* mol : molecules) {
            TrackedMolecule* tracked = find_or_create_tracked(mol, frame_id);
            for (int atom_id : mol->atom_ids) atom_to_molecule[atom_id] = tracked;
        }
        last_processed_frame = frame_id;
        return;
    }

    std::unordered_set<int> changed_atoms;
    for (const auto& [atom_id, new_mol] : current_atom_to_mol) {
        auto old_it = atom_to_molecule.find(atom_id);
        if (old_it != atom_to_molecule.end()) {
            if (old_it->second->structure_id != new_mol->structure_id) {
                changed_atoms.insert(atom_id);
            }
        } else {
            changed_atoms.insert(atom_id);
        }
    }

    if (!changed_atoms.empty()) {
        detect_reactions_from_changes(frame_id, current_atom_to_mol, changed_atoms);
    }

    atom_to_molecule.clear();
    for (Molecule* mol : molecules) {
        TrackedMolecule* tracked = find_or_create_tracked(mol, frame_id);
        for (int atom_id : mol->atom_ids) atom_to_molecule[atom_id] = tracked;
    }

    std::unordered_set<TrackedMolecule*> current_mols;
    for (Molecule* mol : molecules) current_mols.insert(find_or_create_tracked(mol, frame_id));

    for (auto* tracked : active_molecules) {
        if (!current_mols.count(tracked) && tracked->is_active()) {
            tracked->end_frame = frame_id;
        }
    }
    active_molecules = std::move(current_mols);
    last_processed_frame = frame_id;
}

void ReactionTracker::finalize(int last_frame) {
    if (!FLAG_TRACK_REACTIONS) return;
    for (auto* mol : active_molecules) {
        if (mol->is_active()) mol->end_frame = last_frame;
    }
}

void ReactionTracker::save_raw_events(const std::string& filepath) {
    FILE* fp = create_file(filepath);
    if (fp == nullptr) return;

    fmt::print(fp, "event_id,frame,n_reactants,n_products,reactant_ids,product_ids,reactant_formulas,product_formulas,reactant_smiles,product_smiles\n");
    for (const auto& event : events) {
        fmt::print(fp, "{},{},{},{},{},{},{},{},{},{}\n",
                   event.event_id,
                   event.frame_id,
                   event.reactants.size(),
                   event.products.size(),
                   join_ids(event.reactants),
                   join_ids(event.products),
                   csv_escape(join_formulas(event.reactants)),
                   csv_escape(join_formulas(event.products)),
                   csv_escape(join_smiles(event.reactants)),
                   csv_escape(join_smiles(event.products)));
    }
    fclose(fp);
    fmt::print("Saved {} raw reaction events to {}\n", events.size(), filepath);
}

void ReactionTracker::save_raw_event_pairs(const std::string& filepath) {
    FILE* fp = create_file(filepath);
    if (fp == nullptr) return;

    fmt::print(fp, "event_id,frame,source_id,target_id,source_formula,target_formula,source_smiles,target_smiles,tracked_source_id,tracked_target_id,atom_overlap\n");
    int pair_count = 0;
    for (const auto& event : events) {
        for (const auto* reactant : event.reactants) {
            for (const auto* product : event.products) {
                int overlap = 0;
                const auto& smaller = reactant->atom_ids.size() < product->atom_ids.size()
                    ? reactant->atom_ids : product->atom_ids;
                const auto& larger = reactant->atom_ids.size() < product->atom_ids.size()
                    ? product->atom_ids : reactant->atom_ids;
                for (int atom_id : smaller) {
                    if (larger.count(atom_id)) ++overlap;
                }
                if (overlap <= 0) continue;

                fmt::print(fp, "{},{},{},{},{},{},{},{},{},{},{}\n",
                           event.event_id,
                           event.frame_id,
                           reactant->structure_id,
                           product->structure_id,
                           csv_escape(reactant->formula),
                           csv_escape(product->formula),
                           csv_escape(reactant->smiles),
                           csv_escape(product->smiles),
                           reactant->id,
                           product->id,
                           overlap);
                ++pair_count;
            }
        }
    }
    fclose(fp);
    fmt::print("Saved {} raw reaction event atom-overlap pairs to {}\n", pair_count, filepath);
}

void ReactionTracker::save_transfer_flow(const std::string& filepath) {
    struct TransferRecord {
        std::uint64_t source_id = 0;
        std::uint64_t target_id = 0;
        std::string source_formula;
        std::string target_formula;
        std::string source_smiles;
        std::string target_smiles;
        int count = 0;
        int atom_transfer = 0;
    };

    std::map<std::pair<std::uint64_t, std::uint64_t>, TransferRecord> records;
    for (const auto& event : events) {
        for (const auto* reactant : event.reactants) {
            for (const auto* product : event.products) {
                int overlap = 0;
                const auto& smaller = reactant->atom_ids.size() < product->atom_ids.size()
                    ? reactant->atom_ids : product->atom_ids;
                const auto& larger = reactant->atom_ids.size() < product->atom_ids.size()
                    ? product->atom_ids : reactant->atom_ids;
                for (int atom_id : smaller) {
                    if (larger.count(atom_id)) ++overlap;
                }
                if (overlap <= 0) continue;

                auto key = std::make_pair(reactant->structure_id, product->structure_id);
                auto& record = records[key];
                record.source_id = reactant->structure_id;
                record.target_id = product->structure_id;
                record.source_formula = reactant->formula;
                record.target_formula = product->formula;
                record.source_smiles = reactant->smiles;
                record.target_smiles = product->smiles;
                record.count += 1;
                record.atom_transfer += overlap;
            }
        }
    }

    std::vector<TransferRecord> sorted;
    sorted.reserve(records.size());
    for (const auto& [key, record] : records) sorted.push_back(record);
    std::sort(sorted.begin(), sorted.end(), [](const auto& a, const auto& b) {
        if (a.source_id != b.source_id) return a.source_id < b.source_id;
        return a.target_id < b.target_id;
    });

    FILE* fp = create_file(filepath);
    if (fp == nullptr) return;
    fmt::print(fp, "source_id,target_id,source_formula,target_formula,source_smiles,target_smiles,count,atom_transfer,self_loop\n");
    for (const auto& record : sorted) {
        fmt::print(fp, "{},{},{},{},{},{},{},{},{}\n",
                   record.source_id,
                   record.target_id,
                   csv_escape(record.source_formula),
                   csv_escape(record.target_formula),
                   csv_escape(record.source_smiles),
                   csv_escape(record.target_smiles),
                   record.count,
                   record.atom_transfer,
                   record.source_id == record.target_id ? 1 : 0);
    }
    fclose(fp);
    transfer_edge_count = static_cast<int>(sorted.size());
    fmt::print("Saved {} raw transfer-flow edges to {}\n", sorted.size(), filepath);
}

void ReactionTracker::brief_report() const {
    fmt::print("\n=== Reaction Tracker Report ===\n");
    fmt::print("Total molecules tracked: {}\n", all_molecules.size());
    fmt::print("Active molecules: {}\n", active_molecules.size());
    fmt::print("Reaction events detected: {}\n", events.size());

    if (!events.empty()) {
        std::map<std::string, int> counts;
        for (const auto& event : events) {
            std::string key = join_formulas(event.reactants) + " -> " + join_formulas(event.products);
            counts[key] += 1;
        }
        std::vector<std::pair<std::string, int>> sorted(counts.begin(), counts.end());
        std::sort(sorted.begin(), sorted.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        int limit = std::min<int>(10, static_cast<int>(sorted.size()));
        fmt::print("\nTop {} most common reactions:\n", limit);
        for (int i = 0; i < limit; ++i) {
            fmt::print("  {}: {} times\n", sorted[i].first, sorted[i].second);
        }
    }
}
