#include "atom.h"
#include "cell_list.h"
#include "defines.h"
#include "fmt/format.h"
#include "molecule.h"
#include "universe.h"
#include "vec_algorithms.h"

System::System() {}

System::~System() {
    for (auto &atom : atoms) {
        atom->clear();
    }

    for (auto &atom : atoms) {
        delete atom;
    }

    for (auto &bond : bonds) {
        delete bond;
    }

    for (auto &molecule : molecules) {
        delete molecule;
    }

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

void System::finish(std::string &filepath) {
    if (has_boundaries) {
        fmt::print("Atoms: {}, Bonds: {}, Mols: {}, PBC: {:.2f} {:.2f} {:.2f}\n", atoms.size(), bonds.size(),
                   molecules.size(), axis_lengths[0], axis_lengths[1], axis_lengths[2]);
    } else {
        fmt::print("Atoms: {}, Bonds: {}, Mols: {}, PBC: false (cost lots of time)\n", atoms.size(), bonds.size(),
                   molecules.size());
    }
    dump_lammps_data(filepath);
}

// A fallback method for systems without periodic boundaries.
void System::search_neigh_naive(const float &radius, const int &max_neigh) {
    float radius_sq = radius * radius;
    float dist_sq;
    // Naive search for test

    // Define a function pointer type for distance calculation function
    for (auto &curr_atom : atoms) {
        for (auto &other_atom : atoms) {
            if (curr_atom == other_atom) continue;

            dist_sq = distance_sq(curr_atom->coord, other_atom->coord);
            if (dist_sq <= radius) {
                curr_atom->neighs.push_back(other_atom);
            }
            if (curr_atom->neighs.size() >= max_neigh) {
                continue;
            }
        }
    }
}

void System::search_neigh_cell_list(const float &radius, const int &max_neigh) {
    Cell_list cell_list(atoms, radius, axis_lengths, max_neigh);
    for (auto &atom : atoms) {
        cell_list.search_neighbors(atom);
    }
}

void System::search_neigh(const float &radius, const int &max_neigh) {
    if (has_boundaries && axis_lengths.size() == 3 && axis_lengths[0] > 0.0f && axis_lengths[1] > 0.0f &&
        axis_lengths[2] > 0.0f) {
        search_neigh_cell_list(radius, max_neigh);
    } else {
        search_neigh_naive(radius, max_neigh);
    }
}

void System::search_neigh_kdtree(const float &radius, const int &max_neigh) {
    KD_tree kd_tree;

    for (auto &atom : atoms) {
        kd_tree.insert(atom);
    }

    for (auto &atom_data : atoms) {
        if (atom_data->neighs.size() < max_neigh) {
            std::vector<Atom *> neighbors;

            if (has_boundaries && axis_lengths.size() == 3 && axis_lengths[0] > 0.0f && axis_lengths[1] > 0.0f &&
                axis_lengths[2] > 0.0f) {
                kd_tree.find_neighbors(atom_data, neighbors, radius, axis_lengths);
            } else {
                kd_tree.find_neighbors(atom_data, neighbors, radius);
            }

            for (auto &neighbor : neighbors) {
                atom_data->neighs.push_back(neighbor);
                if (atom_data->neighs.size() >= max_neigh) {
                    break;
                }
            }
        }
    }
}

void System::build_bonds_by_radius(const float &rvdw_scale) {
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
        }
    }

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
