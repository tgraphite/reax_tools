#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "cell.h"

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << "\n";
        std::exit(1);
    }
}

bool close(double a, double b, double tol = 1.0e-9) { return std::fabs(a - b) <= tol; }

void test_diagonal() {
    Cell cell;
    cell.set_lattice9({10.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 30.0}, {1.0, 2.0, 3.0},
        {true, true, true});
    require(cell.valid, "diagonal cell valid");
    require(close(cell.a[0], 10.0) && close(cell.b[1], 20.0) && close(cell.c[2], 30.0),
        "diagonal vectors");
    const auto frac = cell.fractional({5.0, 10.0, 15.0});
    require(close(frac[0], 0.5) && close(frac[1], 0.5) && close(frac[2], 0.5), "diagonal fractional");
}

void test_permuted_lattice() {
    Cell cell;
    // This is the OVITO xz-flip form that used to be read as a zero diagonal.
    cell.set_lattice9({0.0, 0.0, 10.0, 0.0, 20.0, 0.0, 30.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
        {true, true, true});
    require(cell.valid, "permuted cell valid");
    require(close(cell.a[2], 10.0) && close(cell.b[1], 20.0) && close(cell.c[0], 30.0),
        "permuted vectors");
    require(cell.is_orthogonal(), "permuted cell is orthogonal");
    const auto frac = cell.fractional({30.0, 20.0, 10.0});
    require(close(frac[0], 1.0) && close(frac[1], 1.0) && close(frac[2], 1.0), "permuted fractional");
}

void test_triclinic() {
    Cell cell;
    const double lx = 10.0, ly = 9.0, lz = 8.0, xy = 3.0, xz = 2.0, yz = 1.0;
    cell.set_lammps_restricted(lx, ly, lz, xy, xz, yz, {0.0, 0.0, 0.0}, {true, true, true});
    require(cell.valid, "triclinic cell valid");
    require(close(cell.a[0], lx) && close(cell.b[0], xy) && close(cell.b[1], ly) &&
                close(cell.c[0], xz) && close(cell.c[1], yz) && close(cell.c[2], lz),
        "triclinic vectors");

    const std::array<double, 3> r{{7.0, 6.0, 5.0}};
    const auto f = cell.fractional(r);
    const auto back = cell.cartesian(f);
    require(close(r[0], back[0]) && close(r[1], back[1]) && close(r[2], back[2]), "fractional roundtrip");

    const std::array<double, 3> delta{{0.6 * lx + 0.6 * xy + 0.6 * xz, 0.6 * ly + 0.6 * yz, 0.6 * lz}};
    const auto minimum = cell.minimum_image(delta);
    const std::array<double, 3> expected{{-0.4 * lx - 0.4 * xy - 0.4 * xz, -0.4 * ly - 0.4 * yz, -0.4 * lz}};
    require(close(minimum[0], expected[0]) && close(minimum[1], expected[1]) && close(minimum[2], expected[2]),
        "minimum image");

    const auto wrapped = cell.wrap_relative({12.0, 11.0, 9.0});
    const auto wf = cell.fractional(wrapped);
    for (int i = 0; i < 3; ++i) {
        require(wf[i] >= -1.0e-9 && wf[i] < 1.0 + 1.0e-9, "wrapped fractional in [0,1)");
    }
}

void test_restricted_conversion() {
    Cell cell;
    cell.set_vectors({2.0, 0.0, 0.0}, {1.0, 3.0, 0.0}, {0.5, 0.7, 4.0}, {0.0, 0.0, 0.0},
        {true, true, true});
    Cell restricted;
    std::array<std::array<double, 3>, 3> rotation{{{{1.0, 0.0, 0.0}}, {{0.0, 1.0, 0.0}}, {{0.0, 0.0, 1.0}}}};
    require(cell.to_restricted(restricted, rotation), "to_restricted succeeds");
    require(restricted.is_restricted(), "restricted form");
    require(close(restricted.a[0], 2.0), "restricted a length");

    const std::array<double, 3> r{{1.2, 0.8, 2.4}};
    std::array<double, 3> rotated{{rotation[0][0] * r[0] + rotation[0][1] * r[1] + rotation[0][2] * r[2],
                                    rotation[1][0] * r[0] + rotation[1][1] * r[1] + rotation[1][2] * r[2],
                                    rotation[2][0] * r[0] + rotation[2][1] * r[1] + rotation[2][2] * r[2]}};
    const double before = Cell::dot(r, r);
    const double after = Cell::dot(rotated, rotated);
    require(close(before, after, 1.0e-9), "restricted rotation preserves distance");

    // Axis flips can produce a left-handed cell; conversion must still work.
    Cell left_handed;
    left_handed.set_lattice9({0.0, 0.0, 10.0, 0.0, 20.0, 0.0, 30.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
        {true, true, true});
    Cell left_restricted;
    std::array<std::array<double, 3>, 3> left_rotation{{{{1.0, 0.0, 0.0}}, {{0.0, 1.0, 0.0}},
                                                         {{0.0, 0.0, 1.0}}}};
    require(left_handed.to_restricted(left_restricted, left_rotation), "left-handed to_restricted");
    require(close(left_restricted.a[0], 10.0) && close(left_restricted.b[1], 20.0) &&
                close(left_restricted.c[2], 30.0),
        "left-handed restricted lengths");
}

}  // namespace

int main() {
    test_diagonal();
    test_permuted_lattice();
    test_triclinic();
    test_restricted_conversion();
    std::cout << "cell_test passed\n";
    return 0;
}
