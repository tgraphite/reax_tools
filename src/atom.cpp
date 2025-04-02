#include "atom.h"

#include <fmt/format.h>

#include <memory>
#include <string>
#include <vector>

// Typically atoms will have less than 12 neighbors in 3.0 A or 3.5 A
// raidus, that is enough for searching topology.
Atom::Atom(int _id, int _type_id, const std::vector<float> &_coord,
           std::string _type_name)
    : id(_id), type_id(_type_id), coord(_coord), type_name(_type_name) {
    neighs.reserve(12);
    bonded_atoms.reserve(8);
}

Atom::~Atom() {
    neighs.clear();
    bonded_atoms.clear();
    bonds.clear();
}

std::string Atom::info() {
    std::string s =
        fmt::format("{:6} {:3} {} {:>6.2f} {:>6.2f} {:>6.2f}\n", id, type_id,
                    type_name, coord[0], coord[1], coord[2]);
    return s;
}

Bond::Bond(std::shared_ptr<Atom> _atom_i, std::shared_ptr<Atom> _atom_j)
    : atom_i(_atom_i), atom_j(_atom_j) {}

Bond::~Bond() {
    atom_i.reset();
    atom_j.reset();
}

std::string Bond::info() {
    std::string s = fmt::format("{}{}-{}{}\n", atom_i->type_name, atom_i->id,
                                atom_j->type_name, atom_j->id);
    return s;
}

Angle::Angle(std::shared_ptr<Atom> _atom_i, std::shared_ptr<Atom> _atom_j,
             std::shared_ptr<Atom> _atom_k)
    : atom_i(_atom_i), atom_j(_atom_j), atom_k(_atom_k) {}

Angle::~Angle() {
    atom_i.reset();
    atom_j.reset();
    atom_k.reset();
}

std::string Angle::info() {
    std::string s = fmt::format("{}{}-{}{}-{}{}\n", atom_i->type_name,
                                atom_i->id, atom_j->type_name, atom_j->id,
                                atom_k->type_name, atom_k->id);
    return s;
}

Dihedral::Dihedral(std::shared_ptr<Atom> _atom_i, std::shared_ptr<Atom> _atom_j,
                   std::shared_ptr<Atom> _atom_k, std::shared_ptr<Atom> _atom_l)
    : atom_i(_atom_i), atom_j(_atom_j), atom_k(_atom_k), atom_l(_atom_l) {}

Dihedral::~Dihedral() {
    atom_i.reset();
    atom_j.reset();
    atom_k.reset();
    atom_k.reset();
}

std::string Dihedral::info() {
    std::string s =
        fmt::format("{}{}-{}{}-{}{}-{}{}\n", atom_i->type_name, atom_i->id,
                    atom_j->type_name, atom_j->id, atom_k->type_name,
                    atom_k->id, atom_l->type_name, atom_l->id);
    return s;
}

// when using (lhs == rhs) here will cause recursive call and blow up stack.
bool operator==(const std::shared_ptr<Atom> &lhs,
                const std::shared_ptr<Atom> &rhs) {
    if (lhs->id == rhs->id)
        return true;
    else
        return false;
}

bool operator==(const std::shared_ptr<Bond> &lhs,
                const std::shared_ptr<Bond> &rhs) {
    // Bond i-j == j-i
    if ((lhs->atom_i == rhs->atom_i) and (lhs->atom_j == rhs->atom_j))
        return true;

    else if ((lhs->atom_i == rhs->atom_j) and (lhs->atom_j == rhs->atom_i))
        return true;

    else
        return false;
}

bool operator==(const std::shared_ptr<Angle> &lhs,
                const std::shared_ptr<Angle> &rhs) {
    // Angle i-j-k == k-j-i
    if (lhs->atom_j != rhs->atom_j) return false;

    if ((lhs->atom_i == rhs->atom_i) and (lhs->atom_k == rhs->atom_k))
        return true;
    else if ((lhs->atom_i == rhs->atom_k) and (lhs->atom_k == rhs->atom_i))
        return true;

    else
        return false;
}

bool operator==(const std::shared_ptr<Dihedral> &lhs,
                const std::shared_ptr<Dihedral> &rhs) {
    // Dihedral i-j-k-l == l-k-j-i
    if ((lhs->atom_i == rhs->atom_i) and (lhs->atom_j == rhs->atom_j) and
        (lhs->atom_k == rhs->atom_k) and (lhs->atom_l == rhs->atom_l))
        return true;

    if ((lhs->atom_i == rhs->atom_l) and (lhs->atom_j == rhs->atom_k) and
        (lhs->atom_k == rhs->atom_j) and (lhs->atom_l == rhs->atom_i))
        return true;

    return false;
}