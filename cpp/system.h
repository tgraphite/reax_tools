#pragma once

#include <fstream>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include "fmt/format.h"
#include "reax_counter.h"
#include "reax_flow.h"
#include "string_tools.h"

class ReaxFlow;
struct Atom;
struct Bond;
struct Molecule;
class System;

struct pair_hash {
    std::size_t operator()(const std::pair<unsigned int, unsigned int>& p) const noexcept {
        return std::hash<unsigned int>()(p.first) ^ (std::hash<unsigned int>()(p.second) << 1);
    }
};

struct Atom {
    // Initialize order. Do not move.
    int id;
    int type_id;
    int max_valence;
    bool saturated = false;
    bool on_ring = false;

    Molecule* belong_molecule = nullptr;

    std::vector<float> coord;
    std::string type_name;
    std::string desc;

    // clear them manually, otherwise atoms will cycle reference each other and
    // the all system will not be released.
    std::unordered_set<Atom*> neighs;
    std::unordered_set<Atom*> bonded_atoms;
    std::unordered_set<Bond*> bonds;

    Atom(int _id, int _type_id, const std::vector<float>& _coord, std::string _type_name)
        : id(_id), type_id(_type_id), coord(_coord), type_name(_type_name) {
        max_valence = ELEMENT_MAX_VALENCIES[type_name];
    };
    ~Atom() {
        neighs.clear();
        bonded_atoms.clear();
        bonds.clear();
        belong_molecule = nullptr;
    };

    inline std::string info() {
        return fmt::format("{:6} {:3} {} {:>6.2f} {:>6.2f} {:>6.2f}\n", id, type_id, type_name, coord[0], coord[1],
            coord[2]);
    };
};

// Make hash of id as atom hash for build unordered_set
namespace std {
    template <>
    struct hash<Atom> {
        std::size_t operator()(const Atom& atom) const noexcept { return std::hash<int>()(atom.id); }
    };
}  // namespace std
inline bool operator==(const Atom& lhs, const Atom& rhs) { return lhs.id == rhs.id; }

struct Bond {
    Atom* atom_i;
    Atom* atom_j;
    int order = 1;

    Bond(Atom* _atom_i, Atom* _atom_j) : atom_i(_atom_i), atom_j(_atom_j) {};
    ~Bond() {};

    inline std::string info() {
        return fmt::format("{}{}-{}{}\n", atom_i->type_name, atom_i->id, atom_j->type_name, atom_j->id);
    }
};

inline bool operator==(const Bond& lhs, const Bond& rhs) {
    return (lhs.atom_i == rhs.atom_i && lhs.atom_j == rhs.atom_j) ||
        (lhs.atom_i == rhs.atom_j && lhs.atom_j == rhs.atom_i);
}

struct Molecule {
    int id;
    unsigned int hash;
    std::string formula;

    std::unordered_set<int> atom_ids;
    std::unordered_set<Atom*> mol_atoms;
    std::unordered_set<Bond*> mol_bonds;
    std::map<std::string, int> types_nums;

    Molecule(int _mol_id) : id(_mol_id) {};
    Molecule(const Molecule& other)
        : id(other.id),
        hash(other.hash),
        atom_ids(other.atom_ids),
        formula(other.formula),
        types_nums(other.types_nums) {
        // Map from old Atom* to new Atom*
        std::map<const Atom*, Atom*> atom_map;

        // Deep copy atoms
        for (const auto& atom_ptr : other.mol_atoms) {
            Atom* new_atom = new Atom(*atom_ptr);
            mol_atoms.insert(new_atom);
            atom_map[atom_ptr] = new_atom;
        }

        // Deep copy bonds, update atom_i/atom_j to new atoms
        for (const auto& bond_ptr : other.mol_bonds) {
            Bond* new_bond = new Bond(*bond_ptr);
            new_bond->atom_i = atom_map[bond_ptr->atom_i];
            new_bond->atom_j = atom_map[bond_ptr->atom_j];
            mol_bonds.insert(new_bond);
        }
    }

    ~Molecule() {
        // Delete all dynamically allocated atoms
        for (auto& atom_ptr : mol_atoms) {
            if (atom_ptr != nullptr) {
                delete atom_ptr;
            }
        }

        // Delete all dynamically allocated bonds
        for (auto& bond_ptr : mol_bonds) {
            if (bond_ptr != nullptr) {
                delete bond_ptr;
            }
        }

        mol_atoms.clear();
        mol_bonds.clear();
    };

    inline std::string info() {
        if (formula.empty()) update_topology();
        return fmt::format("Mol: {} Atoms:{} Bonds:{}\n", formula, mol_atoms.size(), mol_bonds.size());
    };

    inline void update_topology() {
        // Formula
        std::string mol_formula;
        for (auto& mol_atom : mol_atoms) {
            types_nums[mol_atom->type_name]++;
        }

        for (auto& pair : types_nums) {
            mol_formula += fmt::format("{}{}", pair.first, pair.second);
        }
        formula = rename_formula(mol_formula, ELEMENT_DISPLAY_ORDER);

        std::hash<std::string> str_hash;
        hash = str_hash(formula);

        for (auto& atom : mol_atoms) {
            atom->belong_molecule = this;
        }
    };
};

inline bool operator==(const Molecule& lhs, const Molecule& rhs) { return (lhs.hash == rhs.hash); };

class System {
public:
    bool to_destroy = false;
    bool has_boundaries = false;

    bool is_first_frame = false;
    bool is_last_frame = false;

    int total_atoms = 0;
    int total_types = 0;
    int frame_id = 0;

    float neigh_radius = 0.0f;
    unsigned int hash;

    std::vector<float> axis_lengths;
    std::vector<Atom*> atoms;
    std::vector<Molecule*> molecules;
    std::vector<Molecule*> rings;
    std::vector<Bond*> bonds;

    std::unordered_map<int, Atom*> id_atom_map;
    std::unordered_map<int, Molecule*> id_mol_map;

    std::map<std::string, int> type_stoi;
    std::map<int, std::string> type_itos;
    std::map<std::pair<int, int>, float> bond_radius_sq;
    std::map<std::pair<int, int>, float> bond_type_counts;

    std::map<int, int> ring_counts;

    System* prev_sys = nullptr;
    ReaxFlow* reax_flow;
    SpeciesCounter* species_counter;

    Counter<int>* bond_counter;
    Counter<int>* ring_counter;
    Counter<int>* atom_bonded_num_counter;
    Counter<int>* hash_counter;

    System();
    ~System();
    void finish();

    void set_types();
    void load_xyz(std::ifstream& file);
    void load_lammpstrj(std::ifstream& file);

    void search_neigh();
    void search_neigh_naive();
    void search_neigh_cell_list();

    void build_bonds_by_radius();
    void build_molecules();
    void dfs(Atom* atom, std::set<Atom*>& visited, Molecule* curr_mol);

    void process_this();
    void process_counters();
    void process_reax_flow();

    void set_counters(SpeciesCounter* _species_counter, Counter<int>* _bond_counter, Counter<int>* _ring_counter,
        Counter<int>* _atom_bonded_num_counter, Counter<int>* _hash_counter) {
        species_counter = _species_counter;
        bond_counter = _bond_counter;
        ring_counter = _ring_counter;
        atom_bonded_num_counter = _atom_bonded_num_counter;
        hash_counter = _hash_counter;
    };

    void compute_ring_counts();
    void find_rings_from_atom(Atom* current, Atom* start, int depth, std::unordered_set<Atom*>& visited,
        std::unordered_set<std::unordered_set<Atom*>*>& current_rings,
        std::vector<Atom*>& current_path);

    void dump_lammps_data();

    Atom* get_atom_by_id(const int& atom_id);
    Molecule* get_molecule_by_id(const int& mol_id);
};