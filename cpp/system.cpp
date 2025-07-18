#include <assert.h>

#include <algorithm>
#include <mutex>
#include <queue>

#include "cell_list.h"
#include "defines.h"
#include "fmt/format.h"
#include "reax_counter.h"
#include "universe.h"
#include "vec_algorithms.h"

std::mutex reaxcounter_mutex;

// Hash function for std::pair<int, int>
struct pair_hash {
    std::size_t operator()(const std::pair<int, int>& p) const noexcept {
        // Combine the two integers into a single hash value
        // This is a common way to combine two hash values
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
};

System::System() {}

System::~System() {
    // No need to free atoms ptr and bonds ptr, they will be deleted in Molecule destrcutor

    for (const auto& instance : molecules) {
        if (instance != nullptr)
            delete instance;
    }
    molecules.clear();
}

/**
 * @brief Set atom types and initialize type-related mappings and valences.
 * @param type_names Vector of atom type names.
 */
void System::set_types(std::vector<std::string>& type_names) {
    total_types = type_names.size();
    int type_id = 0;
    for (auto& name : type_names) {
        type_id++;
        type_itos[type_id] = name;
        type_stoi[name] = type_id;

        // Initialize ELEMENT_MAX_VALENCIES for each type
        if (ELEMENT_MAX_VALENCIES.find(name) == ELEMENT_MAX_VALENCIES.end()) {
            ELEMENT_MAX_VALENCIES[name] = 4; // Default max valence
        }
    }
}

/**
 * @brief Print system summary for the current frame, including atom, bond, molecule, and ring counts.
 */
void System::finish() {
    fmt::print("\rFrame: {}, Atoms: {}, Bonds: {}, Mols: {}, ", frame_id, atoms.size(), bonds.size(), molecules.size());
    if (has_boundaries) {
        fmt::print("PBC: {:.1f} {:.1f} {:.1f}", axis_lengths[0], axis_lengths[1], axis_lengths[2]);
    } else {
        fmt::print("PBC: false (time++)");
    }

    if (ring_counts.size() > 0) {
        for (auto& ring_count : ring_counts) {
            if (ring_count.second > 0) {
                fmt::print(", R{}={}", ring_count.first, ring_count.second);
            }
        }
    }

    fmt::print("                  "); // If the previous print is too long, clear the line.
}

/**
 * @brief Dump the count of each bond type to a CSV file, combining ij and ji pairs.
 * @param filepath Output file path.
 * @param is_first_frame Indicates if this is the first frame (controls header writing).
 */
void System::dump_bond_count(std::string& filepath, bool& is_first_frame) {
    static FILE* file;

    // consider ij = ji, new_map[ij] = bond_type_counts[ij] + bond_type_counts[ji]
    // std::map<std::pair<int, int>, int> bond_type_counts
    std::map<std::pair<int, int>, int> bond_type_counts_reduced;
    for (const auto& ij_count : bond_type_counts) {
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
        for (auto& ij_count : bond_type_counts_reduced) {
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
    for (auto& ij_count : bond_type_counts_reduced) {
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
 * @brief Dump the count of atoms with different valences to a CSV file.
 * @param filepath Output file path.
 * @param is_first_frame Indicates if this is the first frame (controls header writing).
 */
void System::dump_atom_bonded_num_count(std::string& filepath, bool& is_first_frame) {
    static FILE* file;

    std::map<std::pair<int, int>, int> atom_bonded_num_count; // < atom_type_id, valence >, number

    int type_id;
    std::string type_string;
    int max_valence;

    for (const auto& pair : type_itos) {
        type_id = pair.first;
        type_string = pair.second;

        max_valence = ELEMENT_MAX_VALENCIES[type_string];

        for (int val = 0; val <= max_valence; val++) {
            atom_bonded_num_count[{type_id, val}] = 0;
        }
    }

    // Count
    for (const auto& atom : atoms) {
        atom_bonded_num_count[{atom->type_id, atom->bonded_atoms.size()}]++;
    }

    // Write csv header
    std::string atom_valence_desc;
    if (is_first_frame) {
        file = fopen(filepath.c_str(), "w");
        for (const auto& pair : atom_bonded_num_count) {
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
    for (const auto& pair : atom_bonded_num_count) {
        fmt::print(file, fmt::format("{}", pair.second));
        if (&pair != &*std::prev(atom_bonded_num_count.end())) {
            fmt::print(file, ",");
        }
    }
    fmt::print(file, "\n");
    fclose(file);
}

/**
 * @brief Dump the count of rings of different sizes to a CSV file.
 * @param filepath Output file path.
 * @param is_first_frame Indicates if this is the first frame (controls header writing).
 */
void System::dump_ring_count(std::string& filepath, bool& is_first_frame) {
    static FILE* file;

    if (is_first_frame) {
        file = fopen(filepath.c_str(), "w");
        for (auto& ring_count : ring_counts) {
            fmt::print(file, "size {}", ring_count.first);

            // if ring_count is not the last element.
            if (&ring_count != &*std::prev(ring_counts.end())) {
                fmt::print(file, ",");
            }
        }
        fmt::print(file, "\n");
        fclose(file);
    }

    file = fopen(filepath.c_str(), "a");
    for (auto& ring_count : ring_counts) {
        fmt::print(file, "{}", ring_count.second);

        // if ring_count is not the last element.
        if (&ring_count != &*std::prev(ring_counts.end())) {
            fmt::print(file, ",");
        }
    }
    fmt::print(file, "\n");
    fclose(file);
}

/**
 * @brief Naive neighbor search for systems without periodic boundaries.
 *        Populates each atom's neighbor list based on a cutoff radius.
 */
void System::search_neigh_naive() {
    float radius_sq = neigh_radius * neigh_radius;
    float dist_sq;
    // Naive search for test

    // Define a function pointer type for distance calculation function
    for (auto& curr_atom : atoms) {
        for (auto& other_atom : atoms) {
            if (curr_atom == other_atom)
                continue;

            dist_sq = distance_sq(curr_atom->coord, other_atom->coord);
            if (dist_sq <= radius_sq) {
                curr_atom->neighs.insert(other_atom);
            }
            if (curr_atom->neighs.size() >= max_neigh) {
                continue;
            }
        }
    }
}

/**
 * @brief Use cell list algorithm to search for neighbors within a cutoff radius.
 */
void System::search_neigh_cell_list() {
    Cell_list cell_list(atoms, neigh_radius, axis_lengths, max_neigh);
    for (auto& atom : atoms) {
        cell_list.search_neighbors(atom);
    }
}

/**
 * @brief Selects the appropriate neighbor search method based on boundary conditions.
 */
void System::search_neigh() {
    if (has_boundaries && axis_lengths.size() == 3 && axis_lengths[0] > 0.0f && axis_lengths[1] > 0.0f &&
        axis_lengths[2] > 0.0f) {
        search_neigh_cell_list();
    } else {
        search_neigh_naive();
    }
}

/**
 * @brief Build bonds between atoms based on distance and type-dependent radii.
 *        Ensures valence limits and updates bond/atom connectivity.
 */
void System::build_bonds_by_radius() {
    // prepare types and radius.

    // compute.
    float bond_r;
    float bond_sq;
    float dist_sq;
    float relative_sq;

    Atom* tmp_neigh;
    std::vector<std::pair<Atom*, float>> candidate_id_relative_sq;
    candidate_id_relative_sq.reserve(max_neigh);

    std::pair<int, float> id_dist_sq;
    std::pair<int, int> type_idx;

    for (auto& atom : atoms) {
        // Skip X atoms and atoms that reached max valence
        if (type_itos[atom->type_id] == "X" || atom->bonded_atoms.size() >= atom->max_valence)
            continue;

        for (auto& neigh : atom->neighs) {
            type_idx = {atom->type_id, neigh->type_id};
            bond_sq = bond_radius_sq[type_idx];

            if (bond_sq <= 0)
                continue;

            if (has_boundaries) {
                dist_sq = distance_sq_pbc(atom->coord, neigh->coord, axis_lengths);
            } else {
                dist_sq = distance_sq(atom->coord, neigh->coord);
            }

            if (dist_sq > bond_sq)
                continue;

            relative_sq = dist_sq / bond_sq;

            if (relative_sq < 1.0f) {
                id_dist_sq = {neigh->id, dist_sq};
                candidate_id_relative_sq.emplace_back(std::pair(neigh, relative_sq));
            }
        }

        std::sort(candidate_id_relative_sq.begin(), candidate_id_relative_sq.end(),
                  [](const auto& a, const auto& b) { return a.second < b.second; });

        for (size_t tmp_id = 0; tmp_id < candidate_id_relative_sq.size(); tmp_id++) {
            if (atom->bonded_atoms.size() >= atom->max_valence)
                break;

            tmp_neigh = candidate_id_relative_sq[tmp_id].first;

            if (atom->bonded_atoms.find(tmp_neigh) != atom->bonded_atoms.end())
                continue;

            if (tmp_neigh->bonded_atoms.size() >= tmp_neigh->max_valence)
                continue;

            Bond* bond = new Bond(atom, tmp_neigh);
            bond_type_counts[std::pair(atom->type_id, tmp_neigh->type_id)]++;

            bonds.emplace_back(bond);
            atom->bonds.insert(bond);
            tmp_neigh->bonds.insert(bond);
            atom->bonded_atoms.insert(tmp_neigh);
            tmp_neigh->bonded_atoms.insert(atom);
        }

        candidate_id_relative_sq.clear();
    }
}

/**
 * @brief Build molecules by grouping bonded atoms and maximizing bond orders within valence constraints.
 */
void System::build_molecules() {
    std::set<Atom*> visited;
    int new_mol_id = 0;

    for (auto& atom : atoms) {
        if (visited.find(atom) == visited.end()) // If not visited
        {
            new_mol_id++;
            Molecule* new_mol = new Molecule(new_mol_id);
            dfs(atom, visited, new_mol);
            molecules.push_back(new_mol);
        }
    }

    // Weighing bonds
    // Maximize total weight of edges (bond orders)
    // Constraint: each atom's valence can not exceeded.
    std::map<Atom*, int> atom_valence_residual;
    for (auto& atom : atoms) {
        atom_valence_residual[atom] = atom->max_valence - atom->bonded_atoms.size();
    }

    int current_weight_gain = 0;

    for (auto& molecule : molecules) {
        // Skip single atom molecules and molecules containing X atoms
        if (molecule->mol_atoms.size() == 1)
            continue;

        // Greedy algorithm.
        while (true) {
            current_weight_gain = 0;

            for (auto& bond : molecule->mol_bonds) {
                if (atom_valence_residual[bond->atom_i] > 0 && atom_valence_residual[bond->atom_j] > 0) {
                    current_weight_gain += 2;
                    bond->order++;
                    atom_valence_residual[bond->atom_i]--;
                    atom_valence_residual[bond->atom_j]--;

                    // We do not consider 4+ bonds, they merely exist.
                    if (bond->order >= 3)
                        continue;
                }
            }

            if (current_weight_gain == 0) {
                break;
            }
        }
    }

    atom_valence_residual.clear();

    for (auto& molecule : molecules) {
        molecule->update_topology();
        id_mol_map[molecule->id] = molecule;
    }

    visited.clear();
}

/**
 * @brief Depth-first search to collect all atoms and bonds in a molecule starting from a given atom.
 * @param atom Starting atom for DFS.
 * @param visited Set of already visited atoms.
 * @param curr_mol Molecule object being constructed.
 */
void System::dfs(Atom* atom, std::set<Atom*>& visited, Molecule* curr_mol) {
    visited.insert(atom);
    curr_mol->mol_atoms.insert(atom);
    curr_mol->atom_ids.insert(atom->id);

    for (auto& bond : atom->bonds) {
        curr_mol->mol_bonds.insert(bond);
        // One of atom_i and atom_j is the atom itself.
        if (visited.find(bond->atom_i) == visited.end()) {
            dfs(bond->atom_i, visited, curr_mol);
            continue;
        }
        if (visited.find(bond->atom_j) == visited.end()) {
            dfs(bond->atom_j, visited, curr_mol);
        }
    }
}

/**
 * @brief Dump the current system as a LAMMPS data file for testing.
 * @param filepath Output file path.
 * @param mark_ring_atoms Use mol id entry to mark if an atom on ring, 0 = no , 1 = yes.
 */
void System::dump_lammps_data(std::string& filepath, int frame_step, bool mark_ring_atoms) {
    if (frame_id % frame_step != 0)
        return;

    FILE* file = fopen(filepath.c_str(), "w");
    if (!file) {
        std::cerr << "Failed to open file: " << std::endl;
    }

    fmt::print(file, "# LAMMPS data file written by reax_tools\n\n");
    fmt::print(file, "{} atoms\n", atoms.size());
    fmt::print(file, "{} bonds\n", bonds.size());
    fmt::print(file, "{} atom types\n", total_types);
    fmt::print(file, "{} bond types\n\n", 1);

    if (has_boundaries) {
        fmt::print(file, "{} {} xlo xhi\n", 0, axis_lengths[0]);
        fmt::print(file, "{} {} ylo yhi\n", 0, axis_lengths[1]);
        fmt::print(file, "{} {} zlo zhi\n", 0, axis_lengths[2]);
    }

    fmt::print(file, "\nAtoms # full\n\n");

    if (mark_ring_atoms) {
        for (auto& atom : atoms) {
            fmt::print(file, "{} {} {} 0.0 {:>.3f} {:>.3f} {:>.3f} 0 0 0\n", atom->id, atom->on_ring ? 1 : 0,
                       atom->type_id, atom->coord[0], atom->coord[1], atom->coord[2]);
        }

    } else {
        for (auto& mol : molecules) {
            for (auto& atom : mol->mol_atoms) {
                fmt::print(file, "{} {} {} 0.0 {:>.3f} {:>.3f} {:>.3f} 0 0 0\n", atom->id, mol->id, atom->type_id,
                           atom->coord[0], atom->coord[1], atom->coord[2]);
            }
        }
    }

    fmt::print(file, "\nBonds\n\n");

    for (size_t i = 0; i < bonds.size(); i++) {
        fmt::print(file, "{} 1 {} {}\n", i + 1, bonds[i]->atom_i->id, bonds[i]->atom_j->id);
    }

    fmt::print(file, "\n");

    fclose(file);
}

/**
 * @brief Process the current system: neighbor search, bond building, molecule construction, and ring counting.
 */
void System::process_this() {
    // define bond radius for each type.

    std::map<int, float> atomic_radius;
    float bond_radius = 0.0f;

    for (auto& pair : type_stoi) {
        std::string typ_s = pair.first;
        int typ_i = pair.second;

        // pair.first: string type, pair.second: int type
        if (ELEMENT_ATOMIC_RADII.find(typ_s) != ELEMENT_ATOMIC_RADII.end())
            atomic_radius[typ_i] = ELEMENT_ATOMIC_RADII.at(typ_s);
        else
            atomic_radius[typ_i] = ELEMENT_ATOMIC_RADII.at("X");
    }

    for (int type_i = 1; type_i <= total_types; type_i++) {
        for (int type_j = 1; type_j <= total_types; type_j++) {
            std::pair<int, int> pair_ij = {type_i, type_j};
            if (atomic_radius[type_i] == 0.0f || atomic_radius[type_j] == 0.0f) {
                bond_radius = 0.0f;
            } else {
                bond_radius = 0.5f * (atomic_radius[type_i] + atomic_radius[type_j]) * rvdw_scale;
            }
            if (bond_radius > neigh_radius) {
                neigh_radius = bond_radius;
            }
            bond_radius_sq[pair_ij] = bond_radius * bond_radius;
            bond_type_counts[pair_ij] = 0;
        }
    }

    // neigh_radius = 2.5 * rvdw_scale;

    search_neigh();
    build_bonds_by_radius();
    build_molecules();
    compute_ring_counts();
}

/**
 * @brief Process the current system frame for ReaxFF analysis, including formula extraction and reaction tracking.
 *        Thread-safe import of frame formulas and molecule comparison with previous frame.
 */
void System::process_counters() {
    std::vector<std::string> frame_formulas(molecules.size());
    std::map<std::string, int> bond_count_map;
    std::map<std::string, int> ring_count_map;
    std::map<std::string, int> atom_bonded_num_count_map;

    // Formulas, aka reax_species.
    for (size_t i = 0; i < molecules.size(); i++) {
        frame_formulas[i] = molecules[i]->formula;
    }

    // Bonds
    for (const auto& ij_count : bond_type_counts) {
        int i = ij_count.first.first;
        int j = ij_count.first.second;
        int count = ij_count.second;

        std::pair<int, int> sorted_pair = std::make_pair(std::max(i, j), std::min(i, j));
        std::string key = fmt::format("{}-{}", type_itos[sorted_pair.first], type_itos[sorted_pair.second]);

        if (bond_count_map.find(key) == bond_count_map.end()) {
            bond_count_map[key] = count;
        } else {
            bond_count_map[key] += count;
        }
    }

    // Atom bonded nums
    std::map<std::pair<int, int>, int> tmp_atom_bonded_num_count;
    for (const auto& [type_id, type_string] : type_itos) {
        int max_valence = ELEMENT_MAX_VALENCIES[type_string];
        for (size_t valence = 0; valence <= max_valence; valence++) {
            tmp_atom_bonded_num_count[{type_id, valence}] = 0;
        }
    }
    for (const auto& atom : atoms) {
        tmp_atom_bonded_num_count[{atom->type_id, atom->bonded_atoms.size()}]++;
    }

    for (const auto& [id_num, count] : tmp_atom_bonded_num_count) {
        std::string key = fmt::format("{}-{}", type_itos[id_num.first], id_num.second);
        atom_bonded_num_count_map[key] = count;
    }

    // Rings
    for (const auto& [ring_size, count] : ring_counts) {
        std::string key = fmt::format("ring size {}", ring_size);
        ring_count_map[key] = count;
    }

    // Lock the shared reax_counter import operations.
    {
        std::lock_guard<std::mutex> lock(reaxcounter_mutex);
        reax_counter->import_frame_formulas(frame_id, frame_formulas);
        bond_counter->import_frame_values(frame_id, bond_count_map);
        atom_bonded_num_counter->import_frame_values(frame_id, atom_bonded_num_count_map);
        ring_counter->import_frame_values(frame_id, ring_count_map);
    }
}

inline Atom* System::get_atom_by_id(const int& atom_id) {
    // Use id_atom_map to find the atom by its id. Return nullptr if not found.
    auto it = id_atom_map.find(atom_id);
    if (it != id_atom_map.end()) {
        return it->second;
    } else {
        return nullptr;
    }
}

inline Molecule* System::get_molecule_by_id(const int& mol_id) {
    // Use id_atom_map to find the atom by its id. Return nullptr if not found.
    auto it = id_mol_map.find(mol_id);
    if (it != id_mol_map.end()) {
        return it->second;
    } else {
        return nullptr;
    }
}

void System::process_reax_flow() {
    // These computations are safe without lock.
    if (prev_sys == nullptr) // The first frame.
        return;

    static int has_warned = 5;
    if (prev_sys->atoms.size() != atoms.size() && has_warned > 0) {
        fmt::print("**** WARNING **** : Detected ATOM LOST, reactions may be inaccurate!\n");
        has_warned--;
    } else if (prev_sys->atoms.size() != atoms.size() && has_warned == 0) {
        fmt::print("**** WARNING **** : No more warnings will be printed, you'd better know what you are doing.\n");
        has_warned--;
    }

    std::vector<int> intersection;
    std::unordered_set<std::pair<int, int>, pair_hash> reaction_mol_id_pairs;
    std::pair<int, int> mol_id_pair = {0, 0};
    Atom* prev_atom = nullptr;
    Molecule* prev_mol = nullptr;
    Molecule* curr_mol = nullptr;

    for (const auto& curr_atom : this->atoms) {
        prev_atom = prev_sys->get_atom_by_id(curr_atom->id);
        if (!prev_atom)
            continue;

        prev_mol = prev_atom->belong_molecule;
        curr_mol = curr_atom->belong_molecule;

        if (prev_mol == nullptr || curr_mol == nullptr)
            continue;

        if (prev_mol->hash != curr_mol->hash) {
            mol_id_pair = {prev_mol->id, curr_mol->id};
            if (reaction_mol_id_pairs.find(mol_id_pair) == reaction_mol_id_pairs.end()) {
                reaction_mol_id_pairs.insert(mol_id_pair);
            }
        }
    }

    for (const auto& pair : reaction_mol_id_pairs) {
        intersection.clear();

        prev_mol = prev_sys->get_molecule_by_id(pair.first);
        curr_mol = this->get_molecule_by_id(pair.second);

        std::set_intersection(prev_mol->atom_ids.begin(), prev_mol->atom_ids.end(), curr_mol->atom_ids.begin(),
                              curr_mol->atom_ids.end(), back_inserter(intersection));

        if (intersection.size() > 0) {
            reax_flow->add_reaction(this->frame_id, intersection.size(), prev_mol, curr_mol);
        }
    }
}

/**
 * @brief Analyze molecular reactions and update the reaction flow by comparing molecules between frames.
 *        Records reactions based on atom overlap and contribution metrics.
 */
void System::compute_ring_counts() {
    // Initialize ring counts for sizes 3 to MAX_RING_SIZE
    // MAX_RING_SIZE in defines.h / defines.cpp
    for (int i = MIN_RING_SIZE; i <= MAX_RING_SIZE; i++) {
        ring_counts[i] = 0;
    }

    std::unordered_set<Atom*> visited;
    std::unordered_set<std::unordered_set<Atom*>*> current_rings;
    std::vector<Atom*> current_path;

    // Simple ring detection algorithm.
    for (auto& molecule : molecules) {
        visited.clear();
        current_rings.clear();
        current_path.clear();

        // At least 3 atoms to form rings.
        if (molecule->mol_atoms.size() < 3)
            continue;

        // Search from every atom is mandatory for polycyclics like Naphthalene.
        // Actually this is not very time-consuming.
        for (auto& start_atom : molecule->mol_atoms) {
            find_rings_from_atom(start_atom, start_atom, 0, visited, current_rings, current_path);
        }

        // Set on_ring property for all atoms in rings
        for (auto& ring : current_rings) {
            ring_counts[ring->size()]++;
            // Set on_ring = true for all atoms in this ring
            for (auto& atom : *ring) {
                atom->on_ring = true;
            }
            delete ring; // Only count matters.
        }
    }

    visited.clear();
    current_rings.clear();
    current_path.clear();
}

/**
 * @brief Computes the number of rings of different sizes in the system using DFS-based ring detection.
 *        Updates the ring_counts member variable.
 */
void System::find_rings_from_atom(Atom* current, Atom* start, int depth, std::unordered_set<Atom*>& visited,
                                  std::unordered_set<std::unordered_set<Atom*>*>& current_rings,
                                  std::vector<Atom*>& current_path) {
    // Move forward
    // Add current atom to path
    // depth start from 0 and RING_SIZE start from 1.
    current_path.push_back(current);
    visited.insert(current);
    std::unordered_set<std::unordered_set<Atom*>*> to_erase;
    std::unordered_set<Atom*>* this_ring;

    for (Atom* bonded_atom : current->bonded_atoms) {
        // If ring closed and size = MIN - MAX
        if (bonded_atom == start && depth >= MIN_RING_SIZE - 1 && depth <= MAX_RING_SIZE - 1) {
            this_ring = new std::unordered_set<Atom*>(); // Initialize the pointer
            for (Atom* atom : current_path) {
                this_ring->insert(atom);
            }

            // Determine how to handle this new ring.
            bool to_insert = true;
            if (current_rings.size() > 0) {
                for (auto& other_ring : current_rings) {
                    bool is_subset_of_other = true;
                    bool is_subset_of_this = true;

                    // When this ring has an unique atom, this ring is not a subset of the other one.
                    for (auto& this_atom : *this_ring) {
                        if (other_ring->find(this_atom) == other_ring->end()) {
                            is_subset_of_other = false;
                            break;
                        }
                    }

                    // Vice versa
                    for (auto& other_atom : *other_ring) {
                        if (this_ring->find(other_atom) == this_ring->end()) {
                            is_subset_of_this = false;
                            break;
                        }
                    }

                    // If ring is a proper subset of another, we drop the bigger one.
                    if (is_subset_of_other && !is_subset_of_this) { // The other contains this one
                        to_erase.insert(other_ring);
                        // Don't change to_insert here
                    } else if (is_subset_of_this && !is_subset_of_other) { // This one contains the other
                        to_insert = false;
                        break;                                            // No need to check other rings
                    } else if (is_subset_of_this && is_subset_of_other) { // They are the same
                        to_insert = false;
                        break; // No need to check other rings
                    }
                    // else: None of them is subset, may have intersection, but that is what we want.
                    // to_insert stays true
                }
            }

            for (auto& erase_ring : to_erase) {
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
