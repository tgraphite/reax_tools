#pragma once
#include <map>
#include <memory>
#include <set>
#include <unordered_set>
#include <vector>

#include "atom.h"

struct Molecule {
    int id;
    std::unordered_set<int> atom_ids;

    std::vector<Atom*> mol_atoms;
    std::vector<Bond*> mol_bonds;

    std::string formula;
    std::map<std::string, int> types_nums;

    uint32_t hash;

    Molecule(int _mol_id);
    ~Molecule();

    std::string info();

    void insert(Atom*& atom);
    void insert(Bond*& bond);
    bool has(Atom*& atom);
    bool has(Bond*& bond);

    void update_formula();
    void update_topo();
    void calc_bond_order();

    void dump_mol2(const std::string& filename);
};

// TODO, Important:
// Check the system already have the same molecule, thus no need to generate
// topology again.
bool operator==(const Molecule& lhs, const Molecule& rhs);

// Which means a molecule (sub-molecule, ring) contains another.
bool operator>=(const Molecule& lhs, const Molecule& rhs);
