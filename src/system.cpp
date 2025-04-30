#include <mutex>

#include "atom.h"
#include "cell_list.h"
#include "defines.h"
#include "fmt/format.h"
#include "molecule.h"
#include "reax_species.h"
#include "universe.h"
#include "vec_algorithms.h"

std::mutex reaxspecies_mutex;
std::mutex reaxflow_mutex;

System::System() {}

System::~System() {
    // for (auto &atom : atoms) {
    //     atom->clear();
    // }

    // for (auto &atom : atoms) {
    //     delete atom;
    // }

    // for (auto &bond : bonds) {
    //     delete bond;
    // }

    // for (auto &molecule : molecules) {
    //     delete molecule;
    // }

    atoms.clear();
    bonds.clear();
    molecules.clear();
}

void System::set_types(std::vector<std::string> &type_names) {
    itypes = type_names.size();
    int type_id = 0;
    for (auto &name : type_names) {
        type_id++;
        type_itos[type_id] = name;
        type_stoi[name] = type_id;
    }
}

void System::finish() {
    if (has_boundaries) {
        fmt::print("Atoms: {}, Bonds: {}, Mols: {}, PBC: {:.2f} {:.2f} {:.2f}\n", atoms.size(), bonds.size(),
                   molecules.size(), axis_lengths[0], axis_lengths[1], axis_lengths[2]);
    } else {
        fmt::print("Atoms: {}, Bonds: {}, Mols: {}, PBC: false (cost lots of time)\n", atoms.size(), bonds.size(),
                   molecules.size());
    }
}

void System::dump_bond_count(std::string &filepath, bool &is_first_frame) {
    static FILE *file;

    if (is_first_frame) {
        file = fopen(filepath.c_str(), "w");
        for (auto &ij_count : bond_type_counts) {
            fmt::print(file, "{}-{}", type_itos[ij_count.first.first], type_itos[ij_count.first.second]);

            // if ij_count is not the last element.
            if (&ij_count != &*std::prev(bond_type_counts.end())) {
                fmt::print(file, ",");
            }
        }
        fmt::print(file, "\n");
        fclose(file);
    }

    file = fopen(filepath.c_str(), "a");
    for (auto &ij_count : bond_type_counts) {
        fmt::print(file, "{}", ij_count.second);

        // if ij_count is not the last element.
        if (&ij_count != &*std::prev(bond_type_counts.end())) {
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

    // compute.
    float bond_r;
    float bond_sq;
    float dist_sq;

    for (auto &atom : atoms) {
        for (auto &neigh : atom->neighs) {
            if (neigh == atom) {
                continue;
            }

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

    for (auto &molecule : molecules) {
        molecule->update_formula();
    }
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
        throw std::runtime_error("Failed to open file: " + filepath);
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
}

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
    if (prev_sys == nullptr) return;

    for (const auto &prev_mol : prev_sys->molecules) {
        if (prev_mol == nullptr) continue;

        // Ignore single atom molecule.
        if (prev_mol->atom_ids.size() == 1) continue;

        // Find the most similar molecule in current frame.
        Molecule *best_match = nullptr;
        float best_similarity = 0.0f;

        std::vector<int> intersection;
        std::vector<int> union_set;

        // // Optimization: stop searching once a similarity > 0.5 is found.
        // bool found_match = false;

        for (const auto &curr_mol : this->molecules) {
            if (curr_mol == nullptr) continue;

            // Ignore single atom molecule.
            // if (curr_mol->atom_ids.size() == 1) continue;

            // If the formula is the same, consider it the same molecule, skip.
            if (prev_mol->formula == curr_mol->formula) continue;

            // Quick filter: skip if the size difference is too large.
            // if ((curr_mol->atom_ids.size() / prev_mol->atom_ids.size() >= 2) ||
            //     (prev_mol->atom_ids.size() / curr_mol->atom_ids.size() >= 2))
            //     continue;

            // Calculate similarity: intersection / union.
            intersection.clear();
            union_set.clear();

            std::set_intersection(prev_mol->atom_ids.begin(), prev_mol->atom_ids.end(), curr_mol->atom_ids.begin(),
                                  curr_mol->atom_ids.end(), back_inserter(intersection));

            // Quick check: if intersection is empty, similarity is 0.
            if (intersection.empty()) continue;

            // Quick check: if intersection size equals prev_mol size and
            // curr_mol size, then they are the same molecule, not reaction.
            if (intersection.size() == prev_mol->atom_ids.size() && intersection.size() == curr_mol->atom_ids.size())
                continue;  

            std::set_union(prev_mol->atom_ids.begin(), prev_mol->atom_ids.end(), curr_mol->atom_ids.begin(),
                           curr_mol->atom_ids.end(), back_inserter(union_set));

            float similarity = float(intersection.size()) / float(union_set.size());

            // If similarity exceeds the threshold, record and stop searching.
            if (similarity > 0.5) {
                best_match = curr_mol;
                best_similarity = similarity;
                // found_match = true;
                break;  // If a good match is found, stop searching.
            }

            // If similarity exceeds the threshold, record and stop searching.
            else if (similarity > best_similarity) {
                best_match = curr_mol;
                best_similarity = similarity;
            }
        }

        // If found a match and not the same molecule (lower_limit < similarity < 1.0)
        if (best_match && best_similarity >= 0.2 && best_similarity < 1.0) {
            std::lock_guard<std::mutex> lock(reaxflow_mutex);
            reax_flow->add_reaction(this->frame_id, prev_mol, best_match);
        }
    }
}
