#include "cell_list.h"

#include <cmath>
#include <vector>

#include "fmt/format.h"
#include "system.h"
#include "vec_algorithms.h"

Cell_list::~Cell_list() { cells.clear(); }
// Note that all the atoms have been transformed to {0, lx, 0, ly, 0, lz} before.
Cell_list::Cell_list(std::vector<Atom*> atoms, const float& radius, const std::vector<float> axis_lengths,
                     const int& max_neigh) {
    radius_sq = radius * radius;
    this->max_neigh = max_neigh;
    this->axis_lengths = axis_lengths;

    this->cell_nx = floor(axis_lengths[0] / radius);
    this->cell_ny = floor(axis_lengths[1] / radius);
    this->cell_nz = floor(axis_lengths[2] / radius);
    this->cell_nyz = cell_ny * cell_nz;

    this->cell_lx = axis_lengths[0] / cell_nx;
    this->cell_ly = axis_lengths[1] / cell_ny;
    this->cell_lz = axis_lengths[2] / cell_nz;

    this->inv_cell_lx = 1.0f / cell_lx;
    this->inv_cell_ly = 1.0f / cell_ly;
    this->inv_cell_lz = 1.0f / cell_lz;

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

inline std::vector<int> Cell_list::get_cell_index_xyz(Atom* atom) {
    int cell_idx_x = floor(atom->coord[0] / cell_lx);
    int cell_idx_y = floor(atom->coord[1] / cell_ly);
    int cell_idx_z = floor(atom->coord[2] / cell_lz);

    return {cell_idx_x, cell_idx_y, cell_idx_z};
}

#define WRAP_COORD(val, max) ((val) < 0 ? (max) - 1 : (val) >= (max) ? 0 : (val))
inline std::vector<int> Cell_list::get_neighbor_cell_indices_number(Atom* atom) {
    int cell_idx_x = floor(atom->coord[0] / cell_lx);
    int cell_idx_y = floor(atom->coord[1] / cell_ly);
    int cell_idx_z = floor(atom->coord[2] / cell_lz);

    std::vector<int> neighbor_cell_indices;
    neighbor_cell_indices.reserve(27); // 3x3x3 neighbors

    for (int dx = -1; dx <= 1; ++dx) {
        int x = WRAP_COORD(cell_idx_x + dx, cell_nx);
        for (int dy = -1; dy <= 1; ++dy) {
            int y = WRAP_COORD(cell_idx_y + dy, cell_ny);
            for (int dz = -1; dz <= 1; ++dz) {
                int z = WRAP_COORD(cell_idx_z + dz, cell_nz);
                neighbor_cell_indices.push_back(x * cell_nyz + y * cell_nz + z);
            }
        }
    }
    return neighbor_cell_indices;
}

void Cell_list::search_neighbors(Atom* atom) {
    float distance_sq = 0.0f;
    std::vector<int> neighbor_cell_indices = get_neighbor_cell_indices_number(atom);

    for (int neighbor_cell_index : neighbor_cell_indices) {
        for (Atom* candidate_neighbor : cells[neighbor_cell_index]) {
            if (atom->neighs.size() >= max_neigh) {
                break;
            }

            if (candidate_neighbor != atom) {
                distance_sq = distance_sq_pbc(atom->coord, candidate_neighbor->coord, axis_lengths);
                if (distance_sq < radius_sq) {
                    atom->neighs.insert(candidate_neighbor);
                }
            }
        }
    }
}
