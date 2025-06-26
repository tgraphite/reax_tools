#pragma once

#include <fstream>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include "fmt/format.h"
#include "reax_flow.h"
#include "reax_species.h"
#include "string_tools.h"

class ReaxFlow;
struct Atom;
struct Bond;
struct Molecule;
class System;

struct Atom {
    // Initialize order. Do not move.
    int id;
    int type_id;
    int max_valence;
    bool saturated = false;

    std::vector<float> coord;
    std::string type_name;
    std::string desc;

    // clear them manually, otherwise atoms will cycle reference each other and
    // the all system will not be released.
    std::vector<Atom *> neighs;
    std::vector<Atom *> bonded_atoms;
    std::vector<Bond *> bonds;

    Atom(int _id, int _type_id, const std::vector<float> &_coord, std::string _type_name)
        : id(_id), type_id(_type_id), coord(_coord), type_name(_type_name) {
        max_valence = max_valences[type_name];
    };
    ~Atom() {
        neighs.clear();
        bonded_atoms.clear();
        bonds.clear();
    };

    inline bool contains_neighbor(Atom *atom) {
        if (neighs.size() == 0) {
            return false;
        }

        for (Atom *neigh : neighs) {
            if (neigh == atom) {
                return true;
            }
        }
        return false;
    };

    inline std::string info() {
        return fmt::format("{:6} {:3} {} {:>6.2f} {:>6.2f} {:>6.2f}\n", id, type_id, type_name, coord[0], coord[1],
                           coord[2]);
    };
};

inline bool operator==(const Atom &lhs, const Atom &rhs) { return lhs.id == rhs.id; }

struct Bond {
    Atom *atom_i;
    Atom *atom_j;
    int order = 1;

    Bond(Atom *_atom_i, Atom *_atom_j) : atom_i(_atom_i), atom_j(_atom_j) {};
    ~Bond() {};

    inline std::string info() {
        return fmt::format("{}{}-{}{}\n", atom_i->type_name, atom_i->id, atom_j->type_name, atom_j->id);
    }
};

inline bool operator==(const Bond &lhs, const Bond &rhs) {
    return (lhs.atom_i == rhs.atom_i && lhs.atom_j == rhs.atom_j) ||
           (lhs.atom_i == rhs.atom_j && lhs.atom_j == rhs.atom_i);
}

struct Molecule {
    int id;
    unsigned int hash;
    std::string formula;

    std::unordered_set<int> atom_ids;
    std::unordered_set<Atom *> mol_atoms;
    std::unordered_set<Bond *> mol_bonds;
    std::map<std::string, int> types_nums;

    Molecule(int _mol_id) : id(_mol_id) {};
    Molecule(const Molecule &other)
        : id(other.id), atom_ids(other.atom_ids), formula(other.formula), types_nums(other.types_nums) {
        // Map from old Atom* to new Atom*
        std::map<const Atom *, Atom *> atom_map;

        // Deep copy atoms
        for (const auto &atom_ptr : other.mol_atoms) {
            Atom *new_atom = new Atom(*atom_ptr);
            mol_atoms.insert(new_atom);
            atom_map[atom_ptr] = new_atom;
        }

        // Deep copy bonds, update atom_i/atom_j to new atoms
        for (const auto &bond_ptr : other.mol_bonds) {
            Bond *new_bond = new Bond(*bond_ptr);
            new_bond->atom_i = atom_map[bond_ptr->atom_i];
            new_bond->atom_j = atom_map[bond_ptr->atom_j];
            mol_bonds.insert(new_bond);
        }
    }

    ~Molecule() {
        // Delete all dynamically allocated atoms
        for (auto &atom_ptr : mol_atoms) {
            if (atom_ptr != nullptr) {
                delete atom_ptr;
            }
        }

        // Delete all dynamically allocated bonds
        for (auto &bond_ptr : mol_bonds) {
            if (bond_ptr != nullptr) {
                delete bond_ptr;
            }
        }

        mol_atoms.clear();
        mol_bonds.clear();
    };

    inline std::string info() {
        if (formula.empty()) update_formula();
        return fmt::format("Mol: {} Atoms:{} Bonds:{}\n", formula, mol_atoms.size(), mol_bonds.size());
    };

    inline void update_formula() {
        std::string tmp_formula;

        for (auto &mol_atom : mol_atoms) {
            types_nums[mol_atom->type_name]++;  // If key of map not exist, create
                                                // key-value and value++, If exists,
                                                // value++.
        }

        for (auto &pair : types_nums) {
            tmp_formula += fmt::format("{}{}", pair.first, pair.second);
        }

        formula = rename_formula(tmp_formula);
    };

    inline void update_hash() {
        unsigned int atom_hash = 1;
        unsigned int mol_hash = 1;

        // draft mol hash is build by atom environments, like C(C,H,H) + C(C, C, H) + ...
        // do not depend on the order of atom readed.
        // not absolutely unique, not absolutely represent topology, but fast and very rarely make collision.
        for (const auto &atom : mol_atoms) {
            atom_hash = bigger_prime_numbers[atom->type_id];
            for (const auto &bonded : atom->bonded_atoms) {
                atom_hash *= prime_numbers[bonded->type_id];
            }
            mol_hash += atom_hash;
        }

        hash = mol_hash;
    };
};

inline bool operator==(const Molecule &lhs, const Molecule &rhs) { return (lhs.hash == rhs.hash); };

class System {
   public:
    // total_atoms only useful for loading files.
    bool to_destroy = false;
    bool has_boundaries = false;

    int total_atoms = 0;
    int total_types = 0;
    int frame_id = 0;

    int max_neigh;
    float rvdw_scale;
    unsigned int hash;

    std::vector<float> axis_lengths;
    std::vector<Atom *> atoms;
    std::vector<Molecule *> molecules;
    std::vector<Molecule *> rings;
    std::vector<Bond *> bonds;

    std::map<std::string, int> type_stoi;
    std::map<int, std::string> type_itos;
    std::map<std::pair<int, int>, float> bond_radius_sq;

    std::map<std::pair<int, int>, float> bond_type_counts;
    std::map<int, int> ring_counts;

    System *prev_sys = nullptr;
    ReaxFlow *reax_flow;
    ReaxSpecies *reax_species;

    System();
    ~System();
    void finish();

    void set_types(std::vector<std::string> &type_names);
    void load_xyz(std::ifstream &file);
    void load_lammpstrj(std::ifstream &file);

    void dump_lammps_data(std::string &filepath);
    void dump_bond_count(std::string &filepath, bool &is_first_frame);
    void dump_ring_count(std::string &filepath, bool &is_first_frame);
    void dump_atom_bonded_num_count(std::string &filepath, bool &is_first_frame);

    void search_neigh();
    void search_neigh_naive();
    void search_neigh_cell_list();

    void build_bonds_by_radius();
    void build_molecules();
    void dfs(Atom *atom, std::set<Atom *> &visited, Molecule *curr_mol);

    void process_this();
    void process_reax();

    void set_rvdw_scale(float value) { rvdw_scale = value; }
    void set_max_neigh(int value) { max_neigh = value; }
    void set_reax_flow(ReaxFlow *value) { reax_flow = value; }
    void set_reax_species(ReaxSpecies *value) { reax_species = value; }
    void set_frame_id(int value) { frame_id = value; }
    void set_prev_sys(System *value) { prev_sys = value; }

    void compute_ring_counts();
    void find_rings_from_atom(Atom *current, Atom *start, int depth, std::unordered_set<Atom *> &visited,
                              std::unordered_set<std::unordered_set<Atom *> *> &current_rings,
                              std::vector<Atom *> &current_path);
};