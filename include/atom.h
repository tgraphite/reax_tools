#pragma once
#include <vector>
#include <string>
#include <memory>

struct Atom;
struct Bond;
struct Angle;
struct Dihedral;

struct Atom
{
    // Initialize order. Do not move.
    int id;
    int type_id;
    std::vector<float> coord;
    std::string type_name;

    // clear them manually, otherwise atoms will cycle reference each other and the all system will not be released.
    std::vector<std::shared_ptr<Atom>> neighs;
    std::vector<std::shared_ptr<Atom>> bonded_atoms;
    std::vector<std::shared_ptr<Bond>> bonds;

    std::string desc;

    Atom(int _id, int _type_id, const std::vector<float>& _coord, std::string _type_name);
    ~Atom();

    std::string info();
};

struct Bond
{
    std::shared_ptr<Atom> atom_i;
    std::shared_ptr<Atom> atom_j;

    Bond(std::shared_ptr<Atom> _atom_i, std::shared_ptr<Atom> _atom_j);
    ~Bond();

    std::string info();
};

struct Angle
{
    std::shared_ptr<Atom> atom_i;
    std::shared_ptr<Atom> atom_j;
    std::shared_ptr<Atom> atom_k;

    Angle(std::shared_ptr<Atom> _atom_i, std::shared_ptr<Atom> _atom_j, std::shared_ptr<Atom> _atom_k);
    ~Angle();

    std::string info();
};

struct Dihedral
{
    std::shared_ptr<Atom> atom_i;
    std::shared_ptr<Atom> atom_j;
    std::shared_ptr<Atom> atom_k;
    std::shared_ptr<Atom> atom_l;

    Dihedral(std::shared_ptr<Atom> _atom_i, std::shared_ptr<Atom> _atom_j, std::shared_ptr<Atom> _atom_k, std::shared_ptr<Atom> _atom_l);
    ~Dihedral();

    std::string info();
};

bool operator==(const std::shared_ptr<Atom>& lhs, const std::shared_ptr<Atom>& rhs);
bool operator==(const std::shared_ptr<Bond>& lhs, const std::shared_ptr<Bond>& rhs);
bool operator==(const std::shared_ptr<Angle>& lhs, const std::shared_ptr<Angle>& rhs);
bool operator==(const std::shared_ptr<Dihedral>& lhs, const std::shared_ptr<Dihedral>& rhs);