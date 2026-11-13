#include "cell_list.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "fmt/format.h"
#include "system.h"

Cell_list::~Cell_list() { cells.clear(); }

// Atoms are binned in fractional coordinates.  The number of bins along
// fractional axis i is chosen so that one bin is at least r_cut wide in the
// direction normal to the opposite lattice plane.  With that choice a 3x3x3
// stencil is sufficient even for strongly skewed cells.
Cell_list::Cell_list(const std::vector<Atom*>& atoms, const float& radius, const Cell& in_cell,
                     const int& max_neigh) {
    this->radius = radius;
    radius_sq = radius * radius;
    this->max_neigh = max_neigh;
    cell = in_cell;

    cell_nx = 1;
    cell_ny = 1;
    cell_nz = 1;
    if (radius > 0.0f) {
        for (int i = 0; i < 3; ++i) {
            const double reciprocal_norm = Cell::norm({cell.inv[i][0], cell.inv[i][1], cell.inv[i][2]});
            int n = 1;
            if (reciprocal_norm > 0.0) {
                n = static_cast<int>(std::floor(1.0 / (static_cast<double>(radius) * reciprocal_norm)));
            }
            if (n < 1) n = 1;
            if (i == 0) cell_nx = n;
            if (i == 1) cell_ny = n;
            if (i == 2) cell_nz = n;
        }
    }
    cell_nyz = cell_ny * cell_nz;

    cells.resize(static_cast<std::size_t>(cell_nx) * cell_ny * cell_nz);
    const std::size_t reserve_atoms =
        std::max<std::size_t>(1, (atoms.size() * 2) / std::max<std::size_t>(1, cells.size()));
    for (auto& bin : cells) bin.reserve(reserve_atoms);

    for (auto* atom : atoms) {
        std::array<double, 3> f = cell.fractional_vec(atom->coord);
        for (int i = 0; i < 3; ++i) {
            if (cell.periodic[i]) f[i] -= std::floor(f[i]);
        }
        frac[atom] = f;
        const int idx = get_cell_index(atom);
        if (idx >= 0 && idx < static_cast<int>(cells.size())) cells[idx].push_back(atom);
    }
}

int Cell_list::get_cell_index(Atom* atom) const {
    auto it = frac.find(atom);
    if (it == frac.end()) return -1;
    const auto& f = it->second;
    int ix = static_cast<int>(std::floor(f[0] * cell_nx));
    int iy = static_cast<int>(std::floor(f[1] * cell_ny));
    int iz = static_cast<int>(std::floor(f[2] * cell_nz));
    if (ix < 0) ix = 0;
    if (ix >= cell_nx) ix = cell_nx - 1;
    if (iy < 0) iy = 0;
    if (iy >= cell_ny) iy = cell_ny - 1;
    if (iz < 0) iz = 0;
    if (iz >= cell_nz) iz = cell_nz - 1;
    return ix * cell_nyz + iy * cell_nz + iz;
}

std::vector<int> Cell_list::get_neighbor_cell_indices(Atom* atom) const {
    auto it = frac.find(atom);
    if (it == frac.end()) return {};
    const auto& f = it->second;
    int ix = static_cast<int>(std::floor(f[0] * cell_nx));
    int iy = static_cast<int>(std::floor(f[1] * cell_ny));
    int iz = static_cast<int>(std::floor(f[2] * cell_nz));
    if (ix < 0) ix = 0;
    if (ix >= cell_nx) ix = cell_nx - 1;
    if (iy < 0) iy = 0;
    if (iy >= cell_ny) iy = cell_ny - 1;
    if (iz < 0) iz = 0;
    if (iz >= cell_nz) iz = cell_nz - 1;

    std::array<int, 27> indices{};
    int count = 0;
    for (int dx = -1; dx <= 1; ++dx) {
        const int nx = (ix + dx + cell_nx) % cell_nx;
        for (int dy = -1; dy <= 1; ++dy) {
            const int ny = (iy + dy + cell_ny) % cell_ny;
            for (int dz = -1; dz <= 1; ++dz) {
                const int nz = (iz + dz + cell_nz) % cell_nz;
                indices[count++] = nx * cell_nyz + ny * cell_nz + nz;
            }
        }
    }
    std::sort(indices.begin(), indices.begin() + count);
    const int unique_count = static_cast<int>(std::unique(indices.begin(), indices.begin() + count) - indices.begin());
    return std::vector<int>(indices.begin(), indices.begin() + unique_count);
}

void Cell_list::search_neighbors(Atom* atom) {
    const std::vector<int> neighbor_cell_indices = get_neighbor_cell_indices(atom);
    for (int neighbor_cell_index : neighbor_cell_indices) {
        for (Atom* candidate : cells[neighbor_cell_index]) {
            if (atom->neighs.size() >= static_cast<std::size_t>(max_neigh)) return;
            if (candidate == atom) continue;
            if (cell.distance_sq(atom->coord, candidate->coord) < radius_sq) {
                atom->neighs.insert(candidate);
            }
        }
    }
}
