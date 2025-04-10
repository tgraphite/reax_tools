#include <vec_algorithms.h>

std::vector<float> add_vectors(const std::vector<float>& vec_a, const std::vector<float>& vec_b) {
    std::vector<float> vec(vec_a.size());
    if (!(vec_a.size() == vec_b.size())) {
        std::cerr << "Vector size not equal." << std::endl;
    }

    for (size_t i = 0; i < vec_a.size(); i++) {
        vec[i] = vec_a[i] + vec_b[i];
    }

    return vec;
}

std::vector<float> subtract_vectors(const std::vector<float>& vec_a, const std::vector<float>& vec_b) {
    std::vector<float> vec(vec_a.size());
    if (!(vec_a.size() == vec_b.size())) {
        std::cerr << "Vector size not equal." << std::endl;
    }

    for (size_t i = 0; i < vec_a.size(); i++) {
        vec[i] = vec_a[i] - vec_b[i];
    }

    return vec;
}

std::vector<float> vector_scale(const std::vector<float>& vec_a, const float& k) {
    std::vector<float> vec(vec_a.size());

    for (size_t i = 0; i < vec_a.size(); i++) {
        vec[i] = vec_a[i] * k;
    }

    return vec;
}

std::vector<float> vector_add(const std::vector<float>& vec_a, const float& b) {
    std::vector<float> vec(vec_a.size());

    for (size_t i = 0; i < vec_a.size(); i++) {
        vec[i] = vec_a[i] + b;
    }

    return vec;
}

std::vector<float> operator+(const std::vector<float>& lhs, const std::vector<float>& rhs) {
    return add_vectors(lhs, rhs);
}

std::vector<float> operator-(const std::vector<float>& lhs, const std::vector<float>& rhs) {
    return subtract_vectors(lhs, rhs);
}

std::vector<float> operator*(const std::vector<float>& lhs, const float& rhs) { return vector_scale(lhs, rhs); }

std::vector<float> operator/(const std::vector<float>& lhs, const float& rhs) {
    if (rhs == 0) {
        std::cerr << "Division by zero." << std::endl;
    }

    return vector_scale(lhs, 1 / rhs);
}

float distance_sq(const std::vector<float>& vec_a, const std::vector<float>& vec_b) {
    float dist_sq = 0.0f;

    for (size_t i = 0; i < vec_a.size() && i < vec_b.size(); i++) {
        float dx = vec_a[i] - vec_b[i];
        dist_sq += dx * dx;
    }

    return dist_sq;
}

float distance_sq_pbc(const std::vector<float>& vec_a, const std::vector<float>& vec_b,
                      const std::vector<float>& axis_lengths) {
    float dist_sq = 0.0f;

    for (size_t i = 0; i < vec_a.size() && i < vec_b.size(); i++) {
        float dx = vec_a[i] - vec_b[i];

        if (i < axis_lengths.size() && axis_lengths[i] > 0.0f) {
            dx = dx - axis_lengths[i] * floor(dx / axis_lengths[i] + 0.5f);
        }

        dist_sq += dx * dx;
    }

    return dist_sq;
}
