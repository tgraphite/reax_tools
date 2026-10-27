// Test-only self test for the canonical SMILES engine.
#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "smiles.h"

namespace {

struct Builder {
    std::vector<reax_smiles::Atom> atoms;
    std::vector<reax_smiles::Bond> bonds;

    int add(const std::string& element) {
        int index = static_cast<int>(atoms.size());
        atoms.push_back({index, element, 0});
        return index;
    }
    void bond(int i, int j) { bonds.push_back({i, j}); }
    void hydrogens(int heavy, int count) {
        for (int i = 0; i < count; ++i) bond(heavy, add("H"));
    }
};

reax_smiles::Result build_and_run(Builder& builder) {
    return reax_smiles::canonical_smiles(builder.atoms, builder.bonds);
}

Builder chain(int carbons) {
    Builder b;
    std::vector<int> c;
    for (int i = 0; i < carbons; ++i) c.push_back(b.add("C"));
    for (int i = 0; i + 1 < carbons; ++i) b.bond(c[i], c[i + 1]);
    for (int i = 0; i < carbons; ++i) {
        int h = (i == 0 || i == carbons - 1) ? 3 : 2;
        b.hydrogens(c[i], h);
    }
    return b;
}

Builder ring(int carbons) {
    Builder b;
    std::vector<int> c;
    for (int i = 0; i < carbons; ++i) c.push_back(b.add("C"));
    for (int i = 0; i < carbons; ++i) b.bond(c[i], c[(i + 1) % carbons]);
    for (int i = 0; i < carbons; ++i) b.hydrogens(c[i], 2);
    return b;
}

bool permute_and_check(Builder source, const std::string& expected, std::mt19937& rng) {
    for (int trial = 0; trial < 5; ++trial) {
        std::vector<int> order(source.atoms.size());
        std::iota(order.begin(), order.end(), 0);
        std::shuffle(order.begin(), order.end(), rng);
        std::vector<int> inverse(order.size());
        for (std::size_t i = 0; i < order.size(); ++i) inverse[order[i]] = static_cast<int>(i);

        Builder shuffled;
        for (int old : order) shuffled.atoms.push_back({static_cast<int>(shuffled.atoms.size()),
                                                         source.atoms[old].element,
                                                         source.atoms[old].charge});
        for (const auto& bond : source.bonds) {
            shuffled.bonds.push_back({inverse[bond.i], inverse[bond.j]});
        }
        auto result = build_and_run(shuffled);
        if (result.fallback || result.smiles != expected) return false;
    }
    return true;
}

}  // namespace

int main() {
    std::mt19937 rng(20261027);

    struct Case {
        std::string name;
        Builder builder;
        std::string expected_formula;
    };

    std::vector<Case> cases;

    {
        Builder b;
        int c = b.add("C");
        b.hydrogens(c, 4);
        cases.push_back({"methane", b, "CH4"});
    }
    {
        Builder b;
        int c0 = b.add("C");
        int c1 = b.add("C");
        b.bond(c0, c1);
        b.hydrogens(c0, 3);
        b.hydrogens(c1, 3);
        cases.push_back({"ethane", b, "C2H6"});
    }
    {
        Builder b;
        int c0 = b.add("C");
        int c1 = b.add("C");
        int o = b.add("O");
        b.bond(c0, c1);
        b.bond(c1, o);
        b.hydrogens(c0, 3);
        b.hydrogens(c1, 2);
        b.hydrogens(o, 1);
        cases.push_back({"ethanol", b, "C2H6O"});
    }
    {
        Builder b;
        int c0 = b.add("C");
        int c1 = b.add("C");
        int o = b.add("O");
        b.bond(c0, o);
        b.bond(c1, o);
        b.hydrogens(c0, 3);
        b.hydrogens(c1, 3);
        cases.push_back({"dimethyl_ether", b, "C2H6O"});
    }
    {
        Builder b;
        int n = b.add("N");
        b.hydrogens(n, 3);
        cases.push_back({"ammonia", b, "H3N"});
    }
    {
        Builder b;
        int o = b.add("O");
        b.hydrogens(o, 2);
        cases.push_back({"water", b, "H2O"});
    }
    {
        Builder b;
        int c = b.add("C");
        b.hydrogens(c, 3);
        cases.push_back({"methyl_radical", b, "CH3"});
    }
    {
        Builder b;
        int o = b.add("O");
        b.hydrogens(o, 1);
        cases.push_back({"hydroxyl_radical", b, "HO"});
    }
    {
        Builder b;
        b.bond(b.add("H"), b.add("H"));
        cases.push_back({"h2", b, "H2"});
    }
    {
        Builder b;
        b.add("H");
        cases.push_back({"h_atom", b, "H"});
    }
    {
        Builder b;
        b.bond(b.add("O"), b.add("O"));
        cases.push_back({"o2", b, "O2"});
    }
    {
        Builder b;
        b.bond(b.add("N"), b.add("N"));
        cases.push_back({"n2", b, "N2"});
    }
    {
        Builder b;
        std::vector<int> c;
        for (int i = 0; i < 6; ++i) c.push_back(b.add("C"));
        for (int i = 0; i < 6; ++i) b.bond(c[i], c[(i + 1) % 6]);
        for (int i = 0; i < 6; ++i) b.hydrogens(c[i], 1);
        cases.push_back({"benzene_all_single", b, "C6H6"});
    }
    {
        Builder b;
        std::vector<int> c;
        for (int i = 0; i < 6; ++i) c.push_back(b.add("C"));
        for (int i = 0; i < 6; ++i) b.bond(c[i], c[(i + 1) % 6]);
        for (int i = 0; i < 6; ++i) b.hydrogens(c[i], 2);
        cases.push_back({"cyclohexane", b, "C6H12"});
    }
    {
        // Cubane: 8 carbons, cube edges, one H per carbon.
        Builder b;
        std::vector<int> c;
        for (int i = 0; i < 8; ++i) c.push_back(b.add("C"));
        const int edges[12][2] = {{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},
                                  {6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
        for (auto& e : edges) b.bond(c[e[0]], c[e[1]]);
        for (int i = 0; i < 8; ++i) b.hydrogens(c[i], 1);
        cases.push_back({"cubane", b, "C8H8"});
    }
    {
        // Naphthalene skeleton: two fused six-rings, 10 C, 8 H.
        Builder b;
        std::vector<int> c;
        for (int i = 0; i < 10; ++i) c.push_back(b.add("C"));
        const int edges[11][2] = {{0,1},{1,2},{2,3},{3,4},{4,5},{5,0},
                                  {4,6},{6,7},{7,8},{8,9},{9,3}};
        for (auto& e : edges) b.bond(c[e[0]], c[e[1]]);
        for (int i = 0; i < 10; ++i) {
            int h = (i == 3 || i == 4) ? 0 : 1;
            b.hydrogens(c[i], h);
        }
        cases.push_back({"naphthalene_all_single", b, "C10H8"});
    }
    cases.push_back({"hexane", chain(6), "C6H14"});
    cases.push_back({"cyclododecane", ring(12), "C12H24"});
    // The cutoff counts all atoms, including explicit hydrogens:
    // C66H134 = 200 atoms -> SMILES; C67H136 = 203 atoms -> formula fallback.
    cases.push_back({"chain_66_at_200_atoms", chain(66), "C66H134"});
    cases.push_back({"chain_67_over_200_atoms", chain(67), ""});

    int failures = 0;
    for (auto& item : cases) {
        auto result = build_and_run(item.builder);
        if (item.expected_formula.empty()) {
            if (!result.fallback) {
                std::cout << "FAIL " << item.name << " expected fallback\n";
                ++failures;
            } else {
                std::cout << "OK " << item.name << " FALLBACK "
                          << result.error << "\n";
            }
            continue;
        }
        if (result.fallback || result.smiles.empty()) {
            std::cout << "FAIL " << item.name << " fallback: " << result.error << "\n";
            ++failures;
            continue;
        }
        bool stable = permute_and_check(item.builder, result.smiles, rng);
        std::cout << (stable ? "OK " : "FAIL ") << item.name << " "
                  << result.smiles << (stable ? "" : " (permutation mismatch)")
                  << "\n";
        if (!stable) ++failures;
    }

    std::cout << "failures=" << failures << "\n";
    return failures == 0 ? 0 : 1;
}
