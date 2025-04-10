#pragma once

#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "atom.h"
#include "molecule.h"

class System {
   public:
    // iatoms only useful for loading files.
    int iatoms = 0;
    int itypes = 0;
    int itimestep = 0;
    bool has_boundaries = false;
    std::vector<float> axis_lengths;

    std::vector<Atom *> atoms;
    std::vector<Molecule *> molecules;
    std::vector<Bond *> bonds;

    std::map<std::string, int> type_stoi;
    std::map<int, std::string> type_itos;
    std::map<std::pair<int, int>, float> bond_radius;

    System();
    ~System();

    void set_types(std::vector<std::string> &type_names);
    void load_xyz(std::ifstream &file);
    void load_lammpstrj(std::ifstream &file);

    void finish(std::string &filepath);
    void finish();
    void dump_lammps_data(std::string &filepath);

    void search_neigh(const float &radius, const int &max_neigh);
    void search_neigh_naive(const float &radius, const int &max_neigh);
    void search_neigh_cell_list(const float &radius, const int &max_neigh);
    void search_neigh_kdtree(const float &radius, const int &max_neigh);

    void build_bonds_by_radius(const float &rvdw_scale);
    void build_molecules();
    void dfs(Atom *atom, std::set<Atom *> &visited, Molecule *mol);
};