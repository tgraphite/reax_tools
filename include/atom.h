#pragma once
#include <memory>
#include <string>
#include <vector>

struct Atom;
struct Bond;

struct Atom {
    // Initialize order. Do not move.
    int id;
    int type_id;
    std::vector<float> coord;
    std::string type_name;
    std::string desc;

    // clear them manually, otherwise atoms will cycle reference each other and
    // the all system will not be released.
    std::vector<Atom*> neighs;
    std::vector<Atom*> bonded_atoms;
    std::vector<Bond*> bonds;

    Atom(int _id, int _type_id, const std::vector<float>& _coord, std::string _type_name)
        : id(_id), type_id(_type_id), coord(_coord), type_name(_type_name){};
    ~Atom();

    void clear();

    std::string info();
};

struct Bond {
    Atom* atom_i;
    Atom* atom_j;

    Bond(Atom* _atom_i, Atom* _atom_j) : atom_i(_atom_i), atom_j(_atom_j){};
    ~Bond();

    std::string info();
};

bool operator==(const Atom& lhs, const Atom& rhs);
bool operator==(const Bond& lhs, const Bond& rhs);