#include <algorithm>
#include <mutex>
#include <queue>

#include "atom.h"
#include "cell_list.h"
#include "defines.h"
#include "fmt/format.h"
#include "molecule.h"
#include "reax_species.h"
#include "universe.h"
#include "vec_algorithms.h"

std::mutex reaxspecies_mutex;

System::System() {}

System::~System() {
    // No need to free atoms ptr and bonds ptr, they will be deleted in Molecule destrcutor

    for (const auto &instance : molecules) {
        if (instance != nullptr) delete instance;
    }
    molecules.clear();
}

void System::set_types(std::vector<std::string> &type_names) {
    itypes = type_names.size();
    int type_id = 0;
    for (auto &name : type_names) {
        type_id++;
        type_itos[type_id] = name;
        type_stoi[name] = type_id;

        // Initialize max_valences for each type
        if (max_valences.find(name) == max_valences.end()) {
            max_valences[name] = 4;  // Default max valence
        }
    }
}

void System::finish() {
    fmt::print("\rFrame: {}, Atoms: {}, Bonds: {}, Mols: {}, ", frame_id, atoms.size(), bonds.size(), molecules.size());
    if (has_boundaries) {
        fmt::print("PBC: {:.1f} {:.1f} {:.1f}", axis_lengths[0], axis_lengths[1], axis_lengths[2]);
    } else {
        fmt::print("PBC: false (time++)");
    }

    if (ring_counts.size() > 0) {
        for (auto &ring_count : ring_counts) {
            if (ring_count.second > 0) {
                fmt::print(", R{}={}", ring_count.first, ring_count.second);
            }
        }
    }

    fmt::print("                  ");  // If the previous print is too long, clear the line.
}

void System::dump_bond_count(std::string &filepath, bool &is_first_frame) {
    static FILE *file;

    // consider ij = ji, new_map[ij] = bond_type_counts[ij] + bond_type_counts[ji]
    // std::map<std::pair<int, int>, int> bond_type_counts
    std::map<std::pair<int, int>, int> bond_type_counts_reduced;
    for (const auto &ij_count : bond_type_counts) {
        // New key.
        auto i = ij_count.first.first;
        auto j = ij_count.first.second;
        auto key = std::make_pair(std::max(i, j), std::min(i, j));

        auto count = ij_count.second;

        if (bond_type_counts_reduced.find(key) == bond_type_counts_reduced.end()) {
            bond_type_counts_reduced[key] = count;
        } else {
            bond_type_counts_reduced[key] += count;
        }
    }

    if (is_first_frame) {
        file = fopen(filepath.c_str(), "w");
        for (auto &ij_count : bond_type_counts_reduced) {
            fmt::print(file, "{}-{}", type_itos[ij_count.first.first], type_itos[ij_count.first.second]);

            // if ij_count is not the last element.
            if (&ij_count != &*std::prev(bond_type_counts_reduced.end())) {
                fmt::print(file, ",");
            }
        }
        fmt::print(file, "\n");
        fclose(file);
    }

    file = fopen(filepath.c_str(), "a");
    for (auto &ij_count : bond_type_counts_reduced) {
        fmt::print(file, "{}", ij_count.second);

        // if ij_count is not the last element.
        if (&ij_count != &*std::prev(bond_type_counts_reduced.end())) {
            fmt::print(file, ",");
        }
    }
    fmt::print(file, "\n");
    fclose(file);
}

/**
 * @brief dump count of atom in differenet vanlences, e.g.
 *
 * C(3) - C with 3 bonded atoms (like an sp2 carbon in benzene)
 * O(2) - O with 2 bonded atoms (like H2O)
 */
void System::dump_atom_bonded_num_count(std::string &filepath, bool &is_first_frame) {
    static FILE *file;

    std::map<std::pair<int, int>, int> atom_bonded_num_count;  // < atom_type_id, valence >, number

    int type_id;
    std::string type_string;
    int max_valence;

    for (const auto &pair : type_itos) {
        type_id = pair.first;
        type_string = pair.second;

        max_valence = max_valences[type_string];

        for (int val = 0; val <= max_valence; val++) {
            atom_bonded_num_count[{type_id, val}] = 0;
        }
    }

    // Count
    for (const auto &atom : atoms) {
        atom_bonded_num_count[{atom->type_id, atom->bonded_atoms.size()}]++;
    }

    // Write csv header
    std::string atom_valence_desc;
    if (is_first_frame) {
        file = fopen(filepath.c_str(), "w");
        for (const auto &pair : atom_bonded_num_count) {
            atom_valence_desc = fmt::format("{}({})", type_itos[pair.first.first], pair.first.second);
            fmt::print(file, atom_valence_desc);
            if (&pair != &*std::prev(atom_bonded_num_count.end())) {
                fmt::print(file, ",");
            }
        }
        fmt::print(file, "\n");
        fclose(file);
    }

    // Write data body
    file = fopen(filepath.c_str(), "a");
    for (const auto &pair : atom_bonded_num_count) {
        fmt::print(file, fmt::format("{}", pair.second));
        if (&pair != &*std::prev(atom_bonded_num_count.end())) {
            fmt::print(file, ",");
        }
    }
    fmt::print(file, "\n");
    fclose(file);
}

void System::dump_ring_count(std::string &filepath, bool &is_first_frame) {
    static FILE *file;

    if (is_first_frame) {
        file = fopen(filepath.c_str(), "w");
        for (auto &ring_count : ring_counts) {
            fmt::print(file, "{}-membered ring", ring_count.first);

            // if ring_count is not the last element.
            if (&ring_count != &*std::prev(ring_counts.end())) {
                fmt::print(file, ",");
            }
        }
        fmt::print(file, "\n");
        fclose(file);
    }

    file = fopen(filepath.c_str(), "a");
    for (auto &ring_count : ring_counts) {
        fmt::print(file, "{}", ring_count.second);

        // if ring_count is not the last element.
        if (&ring_count != &*std::prev(ring_counts.end())) {
            fmt::print(file, ",");
        }
    }
    fmt::print(file, "\n");
    fclose(file);
}

// A fallback method for systems without periodic boundaries.
void System::search_neigh_naive() {
    float radius_sq = rvdw_scale * rvdw_scale * 2.5 * 2.5;
    float dist_sq;
    // Naive search for test

    // Define a function pointer type for distance calculation function
    for (auto &curr_atom : atoms) {
        for (auto &other_atom : atoms) {
            if (curr_atom == other_atom) continue;

            dist_sq = distance_sq(curr_atom->coord, other_atom->coord);
            if (dist_sq <= radius_sq) {
                curr_atom->neighs.push_back(other_atom);
            }
            if (curr_atom->neighs.size() >= max_neigh) {
                continue;
            }
        }
    }
}

void System::search_neigh_cell_list() {
    float radius = 2.5 * rvdw_scale;

    Cell_list cell_list(atoms, radius, axis_lengths, max_neigh);
    for (auto &atom : atoms) {
        cell_list.search_neighbors(atom);
    }
}

void System::search_neigh() {
    if (has_boundaries && axis_lengths.size() == 3 && axis_lengths[0] > 0.0f && axis_lengths[1] > 0.0f &&
        axis_lengths[2] > 0.0f) {
        search_neigh_cell_list();
    } else {
        search_neigh_naive();
    }
}

void System::build_bonds_by_radius() {
    // prepare types and radius.
    std::map<int, float> atomic_radius;
    /// When using default atomic radius in constant.h
    for (auto &pair : type_stoi) {
        std::string typ_s = pair.first;
        int typ_i = pair.second;

        // pair.first: string type, pair.second: int type
        if (default_atomic_radius.find(typ_s) != default_atomic_radius.end())
            atomic_radius[typ_i] = default_atomic_radius.at(typ_s);
        else
            atomic_radius[typ_i] = default_atomic_radius.at("X");
    }

    for (int type_i = 1; type_i <= itypes; type_i++) {
        for (int type_j = 1; type_j <= itypes; type_j++) {
            std::pair<int, int> pair_ij = {type_i, type_j};
            if (type_itos[type_i] == "X" || type_itos[type_j] == "X") {
                bond_radius[pair_ij] = 0.0f;
            } else {
                bond_radius[pair_ij] = 0.5f * (atomic_radius[type_i] + atomic_radius[type_j]) * rvdw_scale;
            }
            bond_type_counts[pair_ij] = 0;
        }
    }

    // print bond_type and radius
    // if (frame_id == 1) {
    //     for (auto &pair : bond_type_counts) {
    //         fmt::print("{} {} {}\n", type_itos[pair.first.first], type_itos[pair.first.second],
    //                    bond_radius[pair.first]);
    //     }
    // }

    // compute.
    float bond_r;
    float bond_sq;
    float dist_sq;

    for (auto &atom : atoms) {
        // Skip X atoms and atoms that reached max valence
        if (type_itos[atom->type_id] == "X" || atom->bonded_atoms.size() >= atom->max_valence) continue;

        for (auto &neigh : atom->neighs) {
            // Skip self, X atoms, and atoms that reached max valence
            if (neigh == atom || type_itos[neigh->type_id] == "X") continue;
            if (atom->bonded_atoms.size() >= atom->max_valence || neigh->bonded_atoms.size() >= neigh->max_valence)
                continue;

            if (has_boundaries) {
                dist_sq = distance_sq_pbc(atom->coord, neigh->coord, axis_lengths);
            } else {
                dist_sq = distance_sq(atom->coord, neigh->coord);
            }

            std::pair<int, int> id_ij = {atom->type_id, neigh->type_id};
            bond_r = bond_radius[id_ij];
            bond_sq = bond_r * bond_r;

            if (dist_sq <= bond_sq) {
                Bond *bond = new Bond(atom, neigh);
                bond_type_counts[id_ij]++;

                bonds.push_back(bond);
                atom->bonds.push_back(bond);
                neigh->bonds.push_back(bond);
                atom->bonded_atoms.push_back(neigh);
                neigh->bonded_atoms.push_back(atom);
            }
        }
    }
}

void System::build_molecules() {
    std::set<Atom *> visited;
    int new_mol_id = 0;

    for (auto &atom : atoms) {
        if (visited.find(atom) == visited.end())  // If not visited
        {
            new_mol_id++;
            Molecule *new_mol = new Molecule(new_mol_id);
            dfs(atom, visited, new_mol);
            molecules.push_back(new_mol);
        }
    }

    // Weighing bonds
    // Maximize total weight of edges (bond orders)
    // Constraint: each atom's valence can not exceeded.
    std::map<Atom *, int> atom_valence_residual;
    for (auto &atom : atoms) {
        atom_valence_residual[atom] = atom->max_valence - atom->bonded_atoms.size();
    }

    int current_weight_gain = 0;

    for (auto &molecule : molecules) {
        // Skip single atom molecules and molecules containing X atoms
        if (molecule->mol_atoms.size() == 1) continue;

        // Greedy algorithm.
        while (true) {
            current_weight_gain = 0;

            for (auto &bond : molecule->mol_bonds) {
                if (atom_valence_residual[bond->atom_i] > 0 && atom_valence_residual[bond->atom_j] > 0) {
                    current_weight_gain += 2;
                    bond->order++;
                    atom_valence_residual[bond->atom_i]--;
                    atom_valence_residual[bond->atom_j]--;

                    // We do not consider 4+ bonds, they merely exist.
                    if (bond->order >= 3) continue;
                }
            }

            if (current_weight_gain == 0) {
                break;
            }
        }
    }

    atom_valence_residual.clear();

    for (auto &molecule : molecules) {
        molecule->update_formula();
    }

    visited.clear();
}

void System::dfs(Atom *atom, std::set<Atom *> &visited, Molecule *cur_mol) {
    visited.insert(atom);
    cur_mol->insert(atom);

    for (auto &bond : atom->bonds) {
        cur_mol->insert(bond);
        // One of atom_i and atom_j is the atom itself.
        if (visited.find(bond->atom_i) == visited.end()) {
            dfs(bond->atom_i, visited, cur_mol);
            continue;
        }
        if (visited.find(bond->atom_j) == visited.end()) {
            dfs(bond->atom_j, visited, cur_mol);
        }
    }
}

// Dump lammps data file for test.
void System::dump_lammps_data(std::string &filepath) {
    FILE *file = fopen(filepath.c_str(), "w");
    if (!file) {
        std::cerr << "Failed to open file: " << std::endl;
    }

    fmt::print(file, "# LAMMPS data file written by reax_tools\n\n");
    fmt::print(file, "{} atoms\n", atoms.size());
    fmt::print(file, "{} bonds\n", bonds.size());
    fmt::print(file, "{} atom types\n", itypes);
    fmt::print(file, "{} bond types\n\n", 1);

    if (has_boundaries) {
        fmt::print(file, "{} {} xlo xhi\n", 0, axis_lengths[0]);
        fmt::print(file, "{} {} ylo yhi\n", 0, axis_lengths[1]);
        fmt::print(file, "{} {} zlo zhi\n", 0, axis_lengths[2]);
    }

    fmt::print(file, "\nAtoms # full\n\n");

    for (auto &mol : molecules) {
        for (auto &atom : mol->mol_atoms) {
            // Wrap atom coordinates into the simulation box
            // float x = fmod(atom->coord[0], axis_lengths[0]);
            // float y = fmod(atom->coord[1], axis_lengths[1]);
            // float z = fmod(atom->coord[2], axis_lengths[2]);

            // // Ensure coordinates are positive (in range [0, axis_length])
            // if (x < 0) x += axis_lengths[0];
            // if (y < 0) y += axis_lengths[1];
            // if (z < 0) z += axis_lengths[2];

            // charge is 0 for test, maybe modify later.
            fmt::print(file, "{} {} {} 0.0 {:>.3f} {:>.3f} {:>.3f} 0 0 0\n", atom->id, mol->id, atom->type_id,
                       atom->coord[0], atom->coord[1], atom->coord[2]);
        }
    }

    fmt::print(file, "\nBonds\n\n");

    for (size_t i = 0; i < bonds.size(); i++) {
        fmt::print(file, "{} 1 {} {}\n", i + 1, bonds[i]->atom_i->id, bonds[i]->atom_j->id);
    }

    fmt::print(file, "\n");

    fclose(file);
}

void System::process_this() {
    search_neigh();
    build_bonds_by_radius();
    build_molecules();
    compute_ring_counts();
}

/**
 * @brief Process the current system frame for ReaxFF analysis
 *
 * This function performs several key operations:
 * 1. Extracts molecular formulas from the current frame
 * 2. Imports frame formulas into the ReaxFF species database
 * 3. Compares molecules between current and previous frames to track reactions
 *
 * The function uses thread-safe operations when accessing shared resources
 * and implements a similarity-based matching algorithm to identify unchanged
 * molecules between frames.
 *
 * @note This function assumes that the system has already been processed
 *       through process_this() which builds molecules and computes ring counts
 */
void System::process_reax() {
    std::vector<std::string> frame_formulas(molecules.size());
    for (size_t i = 0; i < molecules.size(); i++) {
        frame_formulas[i] = molecules[i]->formula;
    }

    // Lock the shared reax_species import operations.
    {
        std::lock_guard<std::mutex> lock(reaxspecies_mutex);
        reax_species->import_frame_formulas(frame_id, frame_formulas);
    }

    // These computations are safe without lock.
    if (prev_sys == nullptr)  // The first frame.
        return;

    static int has_warned = 5;
    if (prev_sys->atoms.size() != atoms.size() && has_warned > 0) {
        fmt::print("**** WARNING **** : Detected ATOM LOST, reactions may be inaccurate!\n");
        has_warned--;
    } else if (prev_sys->atoms.size() != atoms.size() && has_warned == 0) {
        fmt::print("**** WARNING **** : No more warnings will be printed, you'd better know what you are doing.\n");
        has_warned--;
    }

    // This double loop compuation is not that time-consuming.
    for (const auto &prev_mol : prev_sys->molecules) {
        if (prev_mol == nullptr) continue;

        // Ignore single atom molecule.
        // if (prev_mol->atom_ids.size() == 1) continue;

        // Find the most similar molecule in current frame.

        bool molecule_unchanged = false;
        Molecule *best_match = nullptr;
        float best_similarity = 0.0f;

        std::vector<int> intersection;
        float contribution = 0.0f;
        // std::vector<int> union_set;

        for (auto &curr_mol : this->molecules) {
            // Ignore single atom molecule.
            // if (curr_mol->atom_ids.size() == 1) continue;
            if (prev_mol == nullptr || curr_mol == nullptr) continue;

            // Molecule operator== is for comparing formulas, maybe topology later, but never atom ids.
            // To make sure if a molecule did not react is that the contribution is 1.0.
            // so if one prev_mol == curr_mol that just mean they are same type, not exact who it is. Just skip.
            if (*prev_mol == *curr_mol) continue;

            // ------------------------------------------------------------
            intersection.clear();

            std::set_intersection(prev_mol->atom_ids.begin(), prev_mol->atom_ids.end(), curr_mol->atom_ids.begin(),
                                  curr_mol->atom_ids.end(), back_inserter(intersection));

            // Quick check: if intersection is empty, similarity is 0.
            if (intersection.empty()) continue;

            // how many % of atoms curr_mol inherits from prev_mol
            // contribution = float(intersection.size()) / float(prev_mol->atom_ids.size());

            // ------------------------------------------------------------

            if (intersection.size() > 0) {
                reax_flow->add_reaction(this->frame_id, intersection.size(), prev_mol, curr_mol);
            }
        }
    }
}

/**
 * @brief Analyzes molecular reactions and updates the reaction flow
 *
 * This function processes molecules in the system to identify and record reactions
 * between consecutive frames. It compares molecules to detect changes and potential
 * reactions based on atom overlap and contribution metrics.
 *
 * The function:
 * - Skips single-atom molecules
 * - Compares molecules between frames using atom ID intersections
 * - Calculates contribution as intersection size / current molecule size
 * - Records reactions when contribution >= 0.5
 *
 * @note This function is part of the reaction tracking system and works in conjunction
 * with the ReaxFlow class to maintain reaction history
 */
void System::compute_ring_counts() {
    // Initialize ring counts for sizes 3 to MAX_RING_SIZE
    // MAX_RING_SIZE in defines.h / defines.cpp
    for (int i = 3; i <= MAX_RING_SIZE; i++) {
        ring_counts[i] = 0;
    }

    std::unordered_set<Atom *> visited;
    std::unordered_set<std::unordered_set<Atom *> *> current_rings;
    std::vector<Atom *> current_path;

    // Simple ring detection algorithm.
    for (auto &molecule : molecules) {
        visited.clear();
        current_rings.clear();
        current_path.clear();

        // At least 3 atoms to form rings.
        if (molecule->mol_atoms.size() < 3) continue;

        // Search from every atom is mandatory for polycyclics like Naphthalene.
        // Actually this is not very time-consuming.
        for (auto &start_atom : molecule->mol_atoms) {
            find_rings_from_atom(start_atom, start_atom, 0, visited, current_rings, current_path);
        }

        for (auto &ring : current_rings) {
            ring_counts[ring->size()]++;
            delete ring;  // Only count matters.
        }
    }

    visited.clear();
    current_rings.clear();
    current_path.clear();
}

/**
 * @brief Computes the number of rings of different sizes in the system
 *
 * This function analyzes all molecules in the system to detect rings of sizes 3 to MAX_RING_SIZE.
 * It uses a depth-first search approach to find rings starting from each atom in each molecule.
 * The results are stored in the ring_counts array where ring_counts[i] represents the number
 * of rings containing i atoms.
 *
 * The function handles:
 * - Multiple rings in a single molecule
 * - Polycyclic systems (e.g., naphthalene)
 * - Proper subset relationships between rings
 *
 * @note This function modifies the ring_counts member variable
 */
void System::find_rings_from_atom(Atom *current, Atom *start, int depth, std::unordered_set<Atom *> &visited,
                                  std::unordered_set<std::unordered_set<Atom *> *> &current_rings,
                                  std::vector<Atom *> &current_path) {
    // Move forward
    // Add current atom to path
    // depth start from 0 and RING_SIZE start from 1.
    current_path.push_back(current);
    visited.insert(current);
    std::unordered_set<std::unordered_set<Atom *> *> to_erase;
    std::unordered_set<Atom *> *this_ring;

    for (Atom *bonded_atom : current->bonded_atoms) {
        // If ring closed and depth = 3-8
        if (bonded_atom == start && depth >= 2 && depth <= MAX_RING_SIZE - 1) {
            this_ring = new std::unordered_set<Atom *>();  // Initialize the pointer
            for (Atom *atom : current_path) {
                this_ring->insert(atom);
            }

            // Determine how to handle this new ring.
            bool to_insert = true;
            if (current_rings.size() > 0) {
                for (auto &other_ring : current_rings) {
                    bool is_subset_of_other = true;
                    bool is_subset_of_this = true;

                    // When this ring has an unique atom, this ring is not a subset of the other one.
                    for (auto &this_atom : *this_ring) {
                        if (other_ring->find(this_atom) == other_ring->end()) {
                            is_subset_of_other = false;
                            break;
                        }
                    }

                    // Vice versa
                    for (auto &other_atom : *other_ring) {
                        if (this_ring->find(other_atom) == this_ring->end()) {
                            is_subset_of_this = false;
                            break;
                        }
                    }

                    // If ring is a proper subset of another, we drop the bigger one.
                    if (is_subset_of_other && !is_subset_of_this) {  // The other contains this one
                        to_erase.insert(other_ring);
                        // Don't change to_insert here
                    } else if (is_subset_of_this && !is_subset_of_other) {  // This one contains the other
                        to_insert = false;
                        break;                                             // No need to check other rings
                    } else if (is_subset_of_this && is_subset_of_other) {  // They are the same
                        to_insert = false;
                        break;  // No need to check other rings
                    }
                    // else: None of them is subset, may have intersection, but that is what we want.
                    // to_insert stays true
                }
            }

            for (auto &erase_ring : to_erase) {
                current_rings.erase(erase_ring);
                delete erase_ring;
            }
            to_erase.clear();

            if (to_insert) {
                current_rings.insert(this_ring);
            } else {
                // Free memory if we're not inserting this ring
                delete this_ring;
            }
        } else if (visited.find(bonded_atom) == visited.end() && depth < MAX_RING_SIZE - 1) {
            // recursive, if not visited
            find_rings_from_atom(bonded_atom, start, depth + 1, visited, current_rings, current_path);
        }
    }

    // Backtrack: remove current atom from path and visited set
    current_path.pop_back();
    visited.erase(current);
}
