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

    std::vector<std::shared_ptr<Atom>> mol_atoms;
    std::vector<std::shared_ptr<Bond>> mol_bonds;
    std::vector<std::shared_ptr<Angle>> mol_angles;
    std::vector<std::shared_ptr<Dihedral>> mol_dihedrals;

    std::string formula;
    std::map<std::string, int> types_nums;

    uint32_t hash;

    Molecule(int _mol_id);
    ~Molecule();

    std::string info();

    void insert(std::shared_ptr<Atom>& atom);
    void insert(std::shared_ptr<Bond>& bond);
    void insert(std::shared_ptr<Angle>& angle);
    void insert(std::shared_ptr<Dihedral>& dihedral);
    bool has(std::shared_ptr<Atom>& atom);
    bool has(std::shared_ptr<Bond>& bond);
    bool has(std::shared_ptr<Angle>& angle);
    bool has(std::shared_ptr<Dihedral>& dihedral);
    void update_formula();
    void update_topo();
};

// TODO, Important:
// Check the system already have the same molecule, thus no need to generate
// topology again.
bool operator==(const Molecule& lhs, const Molecule& rhs);
