#pragma once
#include <cmath>
#include <iostream>
#include <vector>

float dot_product(const std::vector<float> &vec_a, const std::vector<float> &vec_b);
float vector_norm(const std::vector<float> &vec);

std::vector<float> add_vectors(const std::vector<float> &vec_a, const std::vector<float> &vec_b);
std::vector<float> subtract_vectors(const std::vector<float> &vec_a, const std::vector<float> &vec_b);
std::vector<float> vector_scale(const std::vector<float> &vec_a, const float &k);
std::vector<float> vector_add(const std::vector<float> &vec_a, const float &b);
std::vector<float> cross_product(const std::vector<float> &vec_a, const std::vector<float> &vec_b);
std::vector<float> vector_normalize(const std::vector<float> &vec);

std::vector<float> operator+(const std::vector<float> &lhs, const std::vector<float> &rhs);
std::vector<float> operator-(const std::vector<float> &lhs, const std::vector<float> &rhs);

std::vector<float> operator+(const std::vector<float> &lhs, const float &rhs);
std::vector<float> operator-(const std::vector<float> &lhs, const float &rhs);
std::vector<float> operator*(const std::vector<float> &lhs, const float &rhs);
std::vector<float> operator*(const float &lhs, const std::vector<float> &rhs);
std::vector<float> operator/(const std::vector<float> &lhs, const float &rhs);

float distance_sq(const std::vector<float> &vec_a, const std::vector<float> &vec_b,
                  const std::vector<float> &axis_lengths);
float distance_sq_pbc(const std::vector<float> &vec_a, const std::vector<float> &vec_b,
                      const std::vector<float> &axis_lengths);