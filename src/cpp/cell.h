#pragma once

#include <array>
#include <cmath>
#include <vector>

// A general 3D simulation cell.  The three lattice vectors are stored as the
// columns of H, so cartesian = H * fractional.  Atom coordinates in System are
// relative to Cell::origin; this keeps the minimum-image and wrapping maths
// independent of the absolute input origin.
struct Cell {
    std::array<double, 3> a{{0.0, 0.0, 0.0}};
    std::array<double, 3> b{{0.0, 0.0, 0.0}};
    std::array<double, 3> c{{0.0, 0.0, 0.0}};
    std::array<double, 3> origin{{0.0, 0.0, 0.0}};
    std::array<bool, 3> periodic{{false, false, false}};

    bool valid = false;

    // Rows of H^{-1}: fractional_i = dot(inv[i], cartesian).
    std::array<std::array<double, 3>, 3> inv{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};

    static inline double dot(const std::array<double, 3>& u, const std::array<double, 3>& v) {
        return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    }

    static inline std::array<double, 3> cross(
        const std::array<double, 3>& u, const std::array<double, 3>& v) {
        return {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
    }

    static inline double norm(const std::array<double, 3>& u) { return std::sqrt(dot(u, u)); }

    inline std::array<std::array<double, 3>, 3> matrix() const {
        return {{{a[0], b[0], c[0]}, {a[1], b[1], c[1]}, {a[2], b[2], c[2]}}};
    }

    inline bool compute_inverse() {
        double det = dot(a, cross(b, c));
        if (!std::isfinite(det) || std::fabs(det) < 1.0e-12) {
            valid = false;
            inv = {{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
            return false;
        }
        const std::array<double, 3> bc = cross(b, c);
        const std::array<double, 3> ca = cross(c, a);
        const std::array<double, 3> ab = cross(a, b);
        inv[0] = {bc[0] / det, bc[1] / det, bc[2] / det};
        inv[1] = {ca[0] / det, ca[1] / det, ca[2] / det};
        inv[2] = {ab[0] / det, ab[1] / det, ab[2] / det};
        valid = true;
        return true;
    }

    inline void set_vectors(const std::array<double, 3>& va, const std::array<double, 3>& vb,
                            const std::array<double, 3>& vc, const std::array<double, 3>& org,
                            const std::array<bool, 3>& pbc) {
        a = va;
        b = vb;
        c = vc;
        origin = org;
        periodic = pbc;
        compute_inverse();
    }

    // Extended XYZ's 9-value Lattice is Fortran/column-major:
    // ax ay az bx by bz cx cy cz.
    inline void set_lattice9(const std::array<double, 9>& v, const std::array<double, 3>& org,
                             const std::array<bool, 3>& pbc) {
        set_vectors({v[0], v[1], v[2]}, {v[3], v[4], v[5]}, {v[6], v[7], v[8]}, org, pbc);
    }

    inline void set_lattice_diagonal(const std::array<double, 3>& lengths, const std::array<double, 3>& org,
                                     const std::array<bool, 3>& pbc) {
        set_vectors({lengths[0], 0.0, 0.0}, {0.0, lengths[1], 0.0}, {0.0, 0.0, lengths[2]}, org, pbc);
    }

    // LAMMPS restricted triclinic: a=(lx,0,0), b=(xy,ly,0), c=(xz,yz,lz).
    inline void set_lammps_restricted(double lx, double ly, double lz, double xy, double xz, double yz,
                                      const std::array<double, 3>& org, const std::array<bool, 3>& pbc) {
        set_vectors({lx, 0.0, 0.0}, {xy, ly, 0.0}, {xz, yz, lz}, org, pbc);
    }

    inline bool any_periodic() const { return periodic[0] || periodic[1] || periodic[2]; }

    inline std::array<double, 3> fractional(const std::array<double, 3>& r) const {
        return {dot(inv[0], r), dot(inv[1], r), dot(inv[2], r)};
    }

    inline std::array<double, 3> fractional_vec(const std::vector<float>& r) const {
        const std::array<double, 3> arr{{static_cast<double>(r[0]), static_cast<double>(r[1]),
                                          static_cast<double>(r[2])}};
        return fractional(arr);
    }

    inline std::array<double, 3> cartesian(const std::array<double, 3>& f) const {
        return {a[0] * f[0] + b[0] * f[1] + c[0] * f[2], a[1] * f[0] + b[1] * f[1] + c[1] * f[2],
                a[2] * f[0] + b[2] * f[1] + c[2] * f[2]};
    }

    inline std::array<double, 3> wrap_relative(const std::array<double, 3>& r) const {
        if (!valid) return r;
        std::array<double, 3> f = fractional(r);
        for (int i = 0; i < 3; ++i) {
            if (periodic[i]) f[i] -= std::floor(f[i]);
        }
        return cartesian(f);
    }

    inline std::array<double, 3> minimum_image(const std::array<double, 3>& delta) const {
        if (!valid) return delta;
        std::array<double, 3> f = fractional(delta);
        for (int i = 0; i < 3; ++i) {
            if (periodic[i]) f[i] -= std::round(f[i]);
        }
        return cartesian(f);
    }

    inline double distance_sq(const std::array<double, 3>& r1, const std::array<double, 3>& r2) const {
        std::array<double, 3> d{{r1[0] - r2[0], r1[1] - r2[1], r1[2] - r2[2]}};
        if (valid && any_periodic()) d = minimum_image(d);
        return dot(d, d);
    }

    inline double distance_sq(const std::vector<float>& r1, const std::vector<float>& r2) const {
        std::array<double, 3> a1{{static_cast<double>(r1[0]), static_cast<double>(r1[1]),
                                  static_cast<double>(r1[2])}};
        std::array<double, 3> a2{{static_cast<double>(r2[0]), static_cast<double>(r2[1]),
                                  static_cast<double>(r2[2])}};
        return distance_sq(a1, a2);
    }

    inline bool is_orthogonal(double tol = 1.0e-8) const {
        const double na = norm(a);
        const double nb = norm(b);
        const double nc = norm(c);
        if (na < tol || nb < tol || nc < tol) return false;
        return std::fabs(dot(a, b)) <= tol * na * nb && std::fabs(dot(a, c)) <= tol * na * nc &&
               std::fabs(dot(b, c)) <= tol * nb * nc;
    }

    inline bool is_restricted(double tol = 1.0e-8) const {
        return std::fabs(a[1]) <= tol && std::fabs(a[2]) <= tol && std::fabs(b[2]) <= tol;
    }

    inline bool is_restricted_orthogonal(double tol = 1.0e-8) const {
        return is_restricted(tol) && std::fabs(b[0]) <= tol && std::fabs(c[0]) <= tol &&
               std::fabs(c[1]) <= tol;
    }

    // Rotate a right-handed general cell into LAMMPS' restricted lower-triangular
    // form.  rotation maps a coordinate relative to origin into that frame.
    inline bool to_restricted(Cell& out, std::array<std::array<double, 3>, 3>& rotation) const {
        if (!valid) return false;

        // OVITO axis flips can produce a left-handed cell.  The lattice is
        // unchanged if one vector is reversed, so use a right-handed working
        // copy and let the returned rotation absorb the sign change.
        Cell work = *this;
        if (dot(a, cross(b, c)) < 0.0) {
            work.c = {-c[0], -c[1], -c[2]};
            work.compute_inverse();
        }

        const double A = norm(work.a);
        if (A < 1.0e-12) return false;
        const double bx = dot(work.b, work.a) / A;
        const double by = norm(cross(work.a, work.b)) / A;
        if (by < 1.0e-12) return false;
        const double cx = dot(work.c, work.a) / A;
        const double cy = (dot(work.b, work.c) - bx * cx) / by;
        const double cz = dot(cross(work.a, work.b), work.c) / (A * by);
        if (cz < -1.0e-8 || std::fabs(cz) < 1.0e-12) return false;

        Cell restricted;
        restricted.set_vectors({A, 0.0, 0.0}, {bx, by, 0.0}, {cx, cy, cz}, {0.0, 0.0, 0.0}, periodic);

        const auto F = restricted.matrix();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                rotation[i][j] = F[i][0] * work.inv[0][j] + F[i][1] * work.inv[1][j] +
                                 F[i][2] * work.inv[2][j];
            }
        }
        out = restricted;
        return true;
    }
};
