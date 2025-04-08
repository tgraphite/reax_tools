#include "molecule.h"

#include <fmt/format.h>

#include "atom.h"
#include "string_tools.h"

Molecule::Molecule(int _mol_id) : id(_mol_id){};

Molecule::~Molecule() {
    mol_atoms.clear();
    mol_bonds.clear();
}

std::string Molecule::info() {
    if ((mol_atoms.size() > 0) and (formula.empty())) update_formula();

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