#include <cmath>

#include "atom.h"
#include "universe.h"
#include "vec_algorithms.h"

void KD_tree::insert(std::shared_ptr<Atom> &atom_data) { root = insert_rec(root, atom_data, 0); }

KD_tree::node *KD_tree::insert_rec(KD_tree::node *node_ptr, std::shared_ptr<Atom> &atom_data, const int &depth) {
    if (!node_ptr) return new node(atom_data, depth % 3);

    int axis = node_ptr->axis;
    if (atom_data->coord[axis] < node_ptr->atom_data->coord[axis])
        node_ptr->left = insert_rec(node_ptr->left, atom_data, depth + 1);
    else
        node_ptr->right = insert_rec(node_ptr->right, atom_data, depth + 1);

    return node_ptr;
}

void KD_tree::find_neighbors(const std::shared_ptr<Atom> &atom_data, const float &radius,
                             std::vector<std::shared_ptr<Atom>> &neighbors, const std::vector<float> &axis_lengths) {
    find_neighbors_rec(root, atom_data, radius, neighbors, axis_lengths);
}

void KD_tree::find_neighbors_rec(KD_tree::node *node_ptr, const std::shared_ptr<Atom> &atom_data, const float &radius,
                                 std::vector<std::shared_ptr<Atom>> &neighbors,
                                 const std::vector<float> &axis_lengths) {
    if (!node_ptr) return;

    float dist_sq = distance_sq_pbc(node_ptr->atom_data->coord, atom_data->coord, axis_lengths);

    if (dist_sq < radius * radius) neighbors.push_back(node_ptr->atom_data);

    int axis = node_ptr->axis;
    float delta = atom_data->coord[axis] - node_ptr->atom_data->coord[axis];

    if (axis < axis_lengths.size() && axis_lengths[axis] > 0.0f) {
        delta = delta - axis_lengths[axis] * floor(delta / axis_lengths[axis] + 0.5f);
    }

    float radius_sq = radius * radius;

    if (delta < 0) {
        find_neighbors_rec(node_ptr->left, atom_data, radius, neighbors, axis_lengths);
        if (delta * delta < radius_sq) find_neighbors_rec(node_ptr->right, atom_data, radius, neighbors, axis_lengths);
    } else {
        find_neighbors_rec(node_ptr->right, atom_data, radius, neighbors, axis_lengths);
        if (delta * delta < radius_sq) find_neighbors_rec(node_ptr->left, atom_data, radius, neighbors, axis_lengths);
    }
}

void KD_tree::find_neighbors(const std::shared_ptr<Atom> &atom_data, const float &radius,
                             std::vector<std::shared_ptr<Atom>> &neighbors) {
    find_neighbors_rec(root, atom_data, radius, neighbors);
}

void KD_tree::find_neighbors_rec(KD_tree::node *node_ptr, const std::shared_ptr<Atom> &atom_data, const float &radius,
                                 std::vector<std::shared_ptr<Atom>> &neighbors) {
    if (!node_ptr) return;

    float dx = node_ptr->atom_data->coord[0] - atom_data->coord[0];
    float dy = node_ptr->atom_data->coord[1] - atom_data->coord[1];
    float dz = node_ptr->atom_data->coord[2] - atom_data->coord[2];
    float dist_sq = dx * dx + dy * dy + dz * dz;

    if (dist_sq < radius * radius) neighbors.push_back(node_ptr->atom_data);

    int axis = node_ptr->axis;
    float delta = atom_data->coord[axis] - node_ptr->atom_data->coord[axis];
    float delta_sq = delta * delta;
    float radius_sq = radius * radius;
    if (delta < 0) {
        find_neighbors_rec(node_ptr->left, atom_data, radius, neighbors);
        if (delta_sq < radius_sq) find_neighbors_rec(node_ptr->right, atom_data, radius, neighbors);
    } else {
        find_neighbors_rec(node_ptr->right, atom_data, radius, neighbors);
        if (delta_sq < radius_sq) find_neighbors_rec(node_ptr->left, atom_data, radius, neighbors);
    }
}

void KD_tree::clear(KD_tree::node *node_ptr) {
    if (node_ptr) {
        clear(node_ptr->left);
        clear(node_ptr->right);
        delete node_ptr;
    }
}

KD_tree::~KD_tree() { clear(root); }