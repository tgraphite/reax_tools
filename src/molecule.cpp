#include "atom.h"
#include "molecule.h"
#include "string_tools.h"

#include <fmt/format.h>

Molecule::Molecule(int _mol_id) : id(_mol_id){};

Molecule::~Molecule() {
  mol_atoms.clear();
  mol_bonds.clear();
  mol_angles.clear();
  mol_dihedrals.clear();
}

std::string Molecule::info() {
  if ((mol_atoms.size() > 0) and (formula.empty()))
    update_formula();

  std::string s = fmt::format(
      "Mol: {} Bonds:{} Angles:{} Dihedrals:{} hash:{}\n", formula,
      mol_bonds.size(), mol_angles.size(), mol_dihedrals.size(), hash);

  return s;
}

void Molecule::insert(std::shared_ptr<Atom> &atom) {
  if (!this->has(atom)) {
    mol_atoms.push_back(atom);
    atom_ids.insert(atom->id);
  }
}

void Molecule::insert(std::shared_ptr<Bond> &bond) {
  if (!this->has(bond))
    mol_bonds.push_back(bond);
}

void Molecule::insert(std::shared_ptr<Angle> &angle) {
  if (!this->has(angle))
    mol_angles.push_back(angle);
}

void Molecule::insert(std::shared_ptr<Dihedral> &dihedral) {
  if (!this->has(dihedral))
    mol_dihedrals.push_back(dihedral);
}

bool Molecule::has(std::shared_ptr<Atom> &atom) {
  for (auto &exist : mol_atoms) {
    if (atom == exist)
      return true;
  }
  return false;
}

bool Molecule::has(std::shared_ptr<Bond> &bond) {
  for (auto &exist : mol_bonds) {
    if (bond == exist)
      return true;
  }
  return false;
}

bool Molecule::has(std::shared_ptr<Angle> &angle) {
  for (auto &exist : mol_angles) {
    if (angle == exist)
      return true;
  }
  return false;
}

bool Molecule::has(std::shared_ptr<Dihedral> &dihedral) {
  for (auto &exist : mol_dihedrals) {
    if (dihedral == exist)
      return true;
  }
  return false;
}

void Molecule::update_formula() {
  std::string tmp_formula;

  for (auto &mol_atom : mol_atoms) {
    types_nums[mol_atom->type_name]++; // If key of map not exist, create
                                       // key-value and value++, If exists,
                                       // value++.
  }

  for (auto &pair : types_nums) {
    tmp_formula += fmt::format("{}{}", pair.first, pair.second);
  }

  formula = rename_formula(tmp_formula);
}

void Molecule::update_topo() {
  // Build angle from bond.
  for (auto &bond : mol_bonds) {
    // n-i-j
    for (auto &neigh : bond->atom_i->bonded_atoms) {
      if (neigh == bond->atom_j)
        continue;
      auto tmp_angle =
          std::make_shared<Angle>(neigh, bond->atom_i, bond->atom_j);
      this->insert(tmp_angle);
    }
    // i-j-n
    for (auto &neigh : bond->atom_j->bonded_atoms) {
      if (neigh == bond->atom_i)
        continue;
      auto tmp_angle =
          std::make_shared<Angle>(bond->atom_i, bond->atom_j, neigh);
      this->insert(tmp_angle);
    }
  }

  // Build dihedral from angle.
  for (auto &angle : mol_angles) {
    // i-j-k-n
    for (auto &neigh : angle->atom_k->bonded_atoms) {
      if ((neigh == angle->atom_j) or (neigh == angle->atom_i))
        continue;
      auto tmp_dihedral = std::make_shared<Dihedral>(
          angle->atom_i, angle->atom_j, angle->atom_k, neigh);
      this->insert(tmp_dihedral);
    }
    // n-i-j-k
    for (auto &neigh : angle->atom_i->bonded_atoms) {
      if ((neigh == angle->atom_j) or (neigh == angle->atom_k))
        continue;
      auto tmp_dihedral = std::make_shared<Dihedral>(
          neigh, angle->atom_i, angle->atom_j, angle->atom_k);
      this->insert(tmp_dihedral);
    }
  }
}