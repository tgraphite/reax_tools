#include "molecule.h"

#include <fmt/format.h>

#include <algorithm>
#include <functional>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "atom.h"
#include "defines.h"
#include "string_tools.h"

Molecule::Molecule(int _mol_id) : id(_mol_id) {};

Molecule::Molecule(const Molecule &other)
    : id(other.id), atom_ids(other.atom_ids), formula(other.formula), types_nums(other.types_nums), hash(other.hash) {
    // Map from old Atom* to new Atom*
    std::map<const Atom *, Atom *> atom_map;

    // Deep copy atoms
    for (const auto &atom_ptr : other.mol_atoms) {
        Atom *new_atom = new Atom(*atom_ptr);
        mol_atoms.push_back(new_atom);
        atom_map[atom_ptr] = new_atom;
    }

    // Deep copy bonds, update atom_i/atom_j to new atoms
    for (const auto &bond_ptr : other.mol_bonds) {
        Bond *new_bond = new Bond(*bond_ptr);
        new_bond->atom_i = atom_map[bond_ptr->atom_i];
        new_bond->atom_j = atom_map[bond_ptr->atom_j];
        mol_bonds.push_back(new_bond);
    }
}

Molecule::~Molecule() {
    // Delete all dynamically allocated atoms
    for (auto &atom_ptr : mol_atoms) {
        if (atom_ptr != nullptr) {
            delete atom_ptr;
            atom_ptr = nullptr;
        }
    }

    // Delete all dynamically allocated bonds
    for (auto &bond_ptr : mol_bonds) {
        if (bond_ptr != nullptr) {
            delete bond_ptr;
            bond_ptr = nullptr;
        }
    }

    mol_atoms.clear();
    mol_bonds.clear();
}

std::string Molecule::info() {
    if ((mol_atoms.size() > 0) && (formula.empty())) update_formula();

    std::string s = fmt::format("Mol: {} Atoms:{} Bonds:{}\n", formula, mol_atoms.size(), mol_bonds.size());

    return s;
}

void Molecule::insert(Atom *&atom) {
    if (!this->has(atom)) {
        mol_atoms.push_back(atom);
        atom_ids.insert(atom->id);
    }
}

void Molecule::insert(Bond *&bond) {
    if (!this->has(bond)) mol_bonds.push_back(bond);
}

bool Molecule::has(Atom *&atom) {
    for (auto &exist : mol_atoms) {
        if (atom == exist) return true;
    }
    return false;
}

bool Molecule::has(Bond *&bond) {
    for (auto &exist : mol_bonds) {
        if (bond == exist) return true;
    }
    return false;
}

void Molecule::update_formula() {
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
}

void Molecule::dump_mol2(const std::string &filepath) {
    FILE *fp = fopen(filepath.c_str(), "w");
    fmt::print(fp, "# Title\n# Created by ReaxTools\n#\n\n#\n#\n\n");
    fmt::print(fp, "@<TRIPOS>MOLECULE\n{}\n", formula);
    fmt::print(fp, "{} {}\n", mol_atoms.size(), mol_bonds.size());
    fmt::print(fp, "SMALL\nNO_CHARGES\n\n\n");
    fmt::print(fp, "@<TRIPOS>ATOM\n");

    std::map<int, int> map_id_innerid;
    for (int i = 0; i < mol_atoms.size(); i++) {
        map_id_innerid[mol_atoms[i]->id] = i + 1;
    }

    for (int i = 0; i < mol_atoms.size(); i++) {
        fmt::print(fp, "{} {} {} {} {} {}\n", i + 1, mol_atoms[i]->type_name, mol_atoms[i]->coord[0],
                   mol_atoms[i]->coord[1], mol_atoms[i]->coord[2], mol_atoms[i]->type_name);
    }

    fmt::print(fp, "@<TRIPOS>BOND\n");
    for (int i = 0; i < mol_bonds.size(); i++) {
        fmt::print(fp, "{} {} {} {}\n", i + 1, map_id_innerid[mol_bonds[i]->atom_i->id],
                   map_id_innerid[mol_bonds[i]->atom_j->id], std::max(mol_bonds[i]->order, 1));
    }
    fmt::print(fp, "\n");
    fclose(fp);
}

// VERY IMPORTANT: change later: compare by topology, not by formula.
bool operator==(const Molecule &lhs, const Molecule &rhs) { return lhs.formula == rhs.formula; }

// VERY IMPORTANT: change later: compare by topology, not by formula.
bool operator!=(const Molecule &lhs, const Molecule &rhs) { return !(lhs == rhs); }
