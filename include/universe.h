#pragma once
#include <algorithm>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include "atom.h"
#include "molecule.h"
#include "reax_flow.h"
#include "reax_species.h"

class KD_tree;
class System;
class ReaxFlow;

class KD_tree {
   public:
    struct node {
        std::shared_ptr<Atom> atom_data;
        node *left;
        node *right;
        int axis;  // 0 for x, 1 for y, 2 for z

        node(std::shared_ptr<Atom> &atom_data, int axis)
            : atom_data(atom_data), left(nullptr), right(nullptr), axis(axis) {}
    };
    ~KD_tree();

    node *root = nullptr;
    void insert(std::shared_ptr<Atom> &atom_data);
    node *insert_rec(node *node_ptr, std::shared_ptr<Atom> &atom_data,
                     const int &depth);
    void find_neighbors(const std::shared_ptr<Atom> &atom_data,
                        const float &radius,
                        std::vector<std::shared_ptr<Atom>> &neighbors);
    void find_neighbors_rec(node *node_ptr,
                            const std::shared_ptr<Atom> &atom_data,
                            const float &radius,
                            std::vector<std::shared_ptr<Atom>> &neighbors);
    void find_neighbors(const std::shared_ptr<Atom> &atom_data,
                        const float &radius,
                        std::vector<std::shared_ptr<Atom>> &neighbors,
                        const std::vector<float> &axis_lengths);
    void find_neighbors_rec(node *node_ptr,
                            const std::shared_ptr<Atom> &atom_data,
                            const float &radius,
                            std::vector<std::shared_ptr<Atom>> &neighbors,
                            const std::vector<float> &axis_lengths);
    void clear(node *node_ptr);
};

class System {
   public:
    // iatoms only useful for loading files.
    int iatoms = 0;
    int itypes = 0;
    int itimestep = 0;

    bool has_boundaries = false;
    std::vector<float> axis_lengths;

    std::vector<std::shared_ptr<Atom>> atoms;
    std::vector<std::shared_ptr<Molecule>> molecules;
    std::vector<std::shared_ptr<Bond>> bonds;
    std::vector<std::shared_ptr<Angle>> angles;
    std::vector<std::shared_ptr<Dihedral>> dihedrals;

    std::map<std::string, int> type_stoi;
    std::map<int, std::string> type_itos;
    std::map<std::pair<int, int>, float> bond_radius;

    System();
    ~System();

    void set_types(const std::vector<std::string> &type_names);
    void load_xyz(std::ifstream &file);
    void load_lammpstrj(std::ifstream &file);

    void summary();
    void basic_info();

    void search_neigh(const float &r, const int &max_neigh);
    void build_bonds_by_radius(const float &rvdw_scale = 1.2f);
    void build_molecules();
    void dfs(std::shared_ptr<Atom> &atom,
             std::set<std::shared_ptr<Atom>> &visited,
             std::shared_ptr<Molecule> &mol);
};

class Universe {
   public:
    std::shared_ptr<System> prev_sys;
    std::shared_ptr<System> curr_sys;

    std::shared_ptr<ReaxSpecies> reax_species = std::make_shared<ReaxSpecies>();
    std::shared_ptr<ReaxFlow> reax_flow = std::make_shared<ReaxFlow>();

    Universe();
    ~Universe();

    float compute_similarity(const std::unordered_set<int> &prev_set,
                             const std::unordered_set<int> &curr_set);
    void update_reax_flow(const std::shared_ptr<System> &prev_sys,
                          const std::shared_ptr<System> &curr_sys,
                          const int &curr_frame);
    void flush();

    void process_traj(const std::string &file_path,
                      const std::vector<std::string> &type_names,
                      const float &rvdw_scale);
};
