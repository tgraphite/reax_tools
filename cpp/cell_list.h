#include <vector>

#include "system.h"

class Cell_list {
   public:
    std::vector<std::vector<Atom*>> cells;
    std::vector<float> axis_lengths;

    float radius, radius_sq;
    float cell_lx, cell_ly, cell_lz;
    float inv_cell_lx, inv_cell_ly, inv_cell_lz;

    int cell_nx, cell_ny, cell_nz, cell_nyz;
    int max_neigh;

    Cell_list(std::vector<Atom*> atoms, const float& radius, const std::vector<float> axis_lengths,
              const float& max_neigh);
    ~Cell_list();

    std::vector<int> get_cell_index_xyz(Atom* atom);
    std::vector<int> get_neighbor_cell_indices_number(Atom* atom);

    void search_neighbors(Atom* atom);
};