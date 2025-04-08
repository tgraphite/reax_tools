#include <vector>

#include "atom.h"

class KD_tree {
   public:
    struct node {
        Atom *atom_data;
        node *left;
        node *right;
        int axis;  // 0 for x, 1 for y, 2 for z

        node(Atom *atom_data, int axis) : atom_data(atom_data), left(nullptr), right(nullptr), axis(axis) {}
    };
    ~KD_tree();

    node *root = nullptr;
    void insert(Atom *atom_data);
    node *insert_rec(node *node_ptr, Atom *atom_data, const int &depth);

    void find_neighbors(Atom *atom_data, std::vector<Atom *> &neighbors, const float &radius);
    void find_neighbors(Atom *atom_data, std::vector<Atom *> &neighbors, const float &radius,
                        const std::vector<float> &axis_lengths);

    void find_neighbors_rec(node *node_ptr, Atom *atom_data, std::vector<Atom *> &neighbors, const float &radius);
    void find_neighbors_rec(node *node_ptr, Atom *atom_data, std::vector<Atom *> &neighbors, const float &radius,
                            const std::vector<float> &axis_lengths);

    void clear(node *node_ptr);
};