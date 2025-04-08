#include "atom.h"

#include <fmt/format.h>

#include <memory>
#include <string>
#include <vector>

// Typically atoms will have less than 12 neighbors in 3.0 A or 3.5 A
// raidus, that is enough for searching topology.

void Atom::clear() {
    neighs.clear();
    bonded_atoms.clear();
    bonds.clear();
}

Atom::~Atom() { clear(); }

std::string Atom::info() {
    std::string s =
        fmt::format("{:6} {:3} {} {:>6.2f} {:>6.2f} {:>6.2f}\n", id, type_id, type_name, coord[0], coord[1], coord[2]);
    return s;
}

Bond::~Bond() {}

std::string Bond::info() {
    std::string s = fmt::format("{}{}-{}{}\n", atom_i->type_name, atom_i->id, atom_j->type_name, atom_j->id);
    return s;
}

bool operator==(const Atom &lhs, const Atom &rhs) { return lhs.id == rhs.id; }

bool operator==(const Bond &lhs, const Bond &rhs) {
    // Bond i-j == j-i
    return (lhs.atom_i == rhs.atom_i && lhs.atom_j == rhs.atom_j) ||
           (lhs.atom_i == rhs.atom_j && lhs.atom_j == rhs.atom_i);
}