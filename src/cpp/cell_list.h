#pragma once

#include <array>
#include <unordered_map>
#include <vector>

#include "cell.h"
#include "system.h"

class Cell_list {
  public:
    std::vector<std::vector<Atom*>> cells;
    std::unordered_map<Atom*, std::array<double, 3>> frac;
    Cell cell;

    float radius, radius_sq;
    int cell_nx, cell_ny, cell_nz, cell_nyz;
    int max_neigh;

    Cell_list(const std::vector<Atom*>& atoms, const float& radius, const Cell& cell, const int& max_neigh);
    ~Cell_list();

    int get_cell_index(Atom* atom) const;
    std::vector<int> get_neighbor_cell_indices(Atom* atom) const;

    void search_neighbors(Atom* atom);
};
