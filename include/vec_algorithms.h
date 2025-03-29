#pragma once
#include <vector>
#include <iostream>

std::vector<float> add_vectors(const std::vector<float>& vec_a, const std::vector<float>& vec_b);
std::vector<float> subtract_vectors(const std::vector<float>& vec_a, const std::vector<float>& vec_b);
std::vector<float> vector_scale(const std::vector<float>& vec_a, const float& k);
std::vector<float> vector_add(const std::vector<float>& vec_a, const float& b);

// TODO: implement the rest of these operators.
std::vector<float> operator+(const std::vector<float>& lhs, const std::vector<float>& rhs);
std::vector<float> operator-(const std::vector<float>& lhs, const std::vector<float>& rhs);

std::vector<float> operator+(const std::vector<float>& lhs, const float& rhs);
std::vector<float> operator-(const std::vector<float>& lhs, const float& rhs);
std::vector<float> operator*(const std::vector<float>& lhs, const float& rhs);
std::vector<float> operator/(const std::vector<float>& lhs, const float& rhs);