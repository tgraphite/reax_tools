#pragma once

#include <fstream>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "atom.h"
#include "molecule.h"
#include "reax_flow.h"
#include "reax_species.h"

class System {
   public:
    // iatoms only useful for loading files.
    int iatoms = 0;
    int itypes = 0;
    int itimestep = 0;
    int frame_id = 0;

    bool has_boundaries = false;
    float rvdw_scale;
    int max_neigh;
    float reaxflow_threshold;

    std::vector<float> axis_lengths;

    std::vector<Atom *> atoms;
    std::vector<Molecule *> molecules;
    std::vector<Molecule *> rings;
    std::vector<Bond *> bonds;

    std::map<std::string, int> type_stoi;
    std::map<int, std::string> type_itos;
    std::map<std::pair<int, int>, float> bond_radius;
    std::map<std::pair<int, int>, float> bond_type_counts;

    std::map<int, int> ring_counts;

    System *prev_sys = nullptr;
    ReaxFlow *reax_flow;
    ReaxSpecies *reax_species;

    System();
    ~System();

    void set_types(std::vector<std::string> &type_names);
    void load_xyz(std::ifstream &file);
    void load_lammpstrj(std::ifstream &file);

    void finish();
    void dump_lammps_data(std::string &filepath);
    void dump_bond_count(std::string &filepath, bool &is_first_frame);

    void search_neigh();
    void search_neigh_naive();
    void search_neigh_cell_list();
    void search_neigh_kdtree();

    void build_bonds_by_radius();
    void build_molecules();
    void dfs(Atom *atom, std::set<Atom *> &visited, Molecule *mol);

    void process_this();
    void process_reax();
    float compute_similarity(std::unordered_set<int> &prev_set, std::unordered_set<int> &curr_set);
    void process_outputs();

    void set_rvdw_scale(float value) { rvdw_scale = value; }
    void set_max_neigh(int value) { max_neigh = value; }
    void set_reax_flow(ReaxFlow *value) { reax_flow = value; }
    void set_reax_species(ReaxSpecies *value) { reax_species = value; }
    void set_frame_id(int value) { frame_id = value; }
    void set_prev_sys(System *value) { prev_sys = value; }
    void set_reaxflow_threshold(float value) { reaxflow_threshold = value; }

    void compute_ring_counts();
    // Helper function to find cycles in molecules using DFS with depth limit
    void find_rings_from_atom(Atom *current, Atom *start, int depth, std::unordered_set<Atom *> &visited,
                              std::unordered_set<std::unordered_set<Atom *> *> &current_rings,
                              std::vector<Atom *> &current_path);
};