#include "cell_list.h"

#include <fmt/core.h>

#include <cmath>
#include <vector>

#include "atom.h"
#include "vec_algorithms.h"

Cell_list::~Cell_list() { cells.clear(); }
// Note that all the atoms have been transformed to {0, lx, 0, ly, 0, lz} before.
Cell_list::Cell_list(std::vector<Atom*> atoms, const float& radius, const std::vector<float> axis_lengths,
                     const float& max_neigh) {
    radius_sq = radius * radius;
    this->max_neigh = max_neigh;
    this->axis_lengths = axis_lengths;

    this->cell_nx = floor(axis_lengths[0] / radius);
    this->cell_ny = floor(axis_lengths[1] / radius);
    this->cell_nz = floor(axis_lengths[2] / radius);

    this->cell_lx = axis_lengths[0] / cell_nx;
    this->cell_ly = axis_lengths[1] / cell_ny;
    this->cell_lz = axis_lengths[2] / cell_nz;

    int cell_idx_x, cell_idx_y, cell_idx_z;

    cells.resize(cell_nx * cell_ny * cell_nz);
    int reserve_atoms = ceil(atoms.size() / cells.size() * 2);
    for (auto& cell : cells) {
        cell.reserve(reserve_atoms);
    }

    for (auto& atom : atoms) {
        std::vector<int> cell_idx_xyz = get_cell_index_xyz(atom);
        int cell_idx = cell_idx_xyz[0] * cell_ny * cell_nz + cell_idx_xyz[1] * cell_nz + cell_idx_xyz[2];

        // Add bounds checking to prevent segfault
        if (cell_idx >= 0 && cell_idx < cells.size()) {
            Atom* copy_ptr = atom;
            cells[cell_idx].push_back(copy_ptr);
        }
    }
}

std::vector<int> Cell_list::get_cell_index_xyz(Atom* atom) {
    int cell_idx_x = floor(atom->coord[0] / cell_lx);
    int cell_idx_y = floor(atom->coord[1] / cell_ly);
    int cell_idx_z = floor(atom->coord[2] / cell_lz);

    return {cell_idx_x, cell_idx_y, cell_idx_z};
}

std::vector<int> Cell_list::get_neighbor_cell_indices_number(Atom* atom) {
    std::vector<int> cell_idx_xyz = get_cell_index_xyz(atom);
    std::vector<int> neighbor_cell_indices;
    for (int move_x = -1; move_x <= 1; move_x++) {
        int new_x = cell_idx_xyz[0] + move_x;
        if (new_x < 0) {
            new_x = cell_nx - 1;
        } else if (new_x >= cell_nx) {
            new_x = 0;
        }
        for (int move_y = -1; move_y <= 1; move_y++) {
            int new_y = cell_idx_xyz[1] + move_y;
            if (new_y < 0) {
                new_y = cell_ny - 1;
            } else if (new_y >= cell_ny) {
                new_y = 0;
            }
            for (int move_z = -1; move_z <= 1; move_z++) {
                int new_z = cell_idx_xyz[2] + move_z;
                if (new_z < 0) {
                    new_z = cell_nz - 1;
                } else if (new_z >= cell_nz) {
                    new_z = 0;
                }
                neighbor_cell_indices.push_back(new_x * cell_ny * cell_nz + new_y * cell_nz + new_z);
            }
        }
    }
    return neighbor_cell_indices;
}

void Cell_list::search_neighbors(Atom* atom) {
    float distance_sq = 10.0;
    std::vector<int> neighbor_cell_indices = get_neighbor_cell_indices_number(atom);
    for (int neighbor_cell_index : neighbor_cell_indices) {
        for (Atom* candidate_neighbor : cells[neighbor_cell_index]) {
            if (atom->neighs.size() >= max_neigh) {
                break;
            }

            if (atom->contains_neighbor(candidate_neighbor) || candidate_neighbor->contains_neighbor(atom)) {
                continue;
            }

            if (candidate_neighbor != atom) {
                distance_sq = distance_sq_pbc(atom->coord, candidate_neighbor->coord, axis_lengths);
                if (distance_sq < radius_sq) {
                    atom->neighs.push_back(candidate_neighbor);
                }
            }
        }
    }
}
