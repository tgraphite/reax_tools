// Test-only CLI: read an explicit-atom graph from stdin and print SMILES.
// Input:
//   <n> <m>
//   <element>            (n lines; charge optional as "<element> <charge>")
//   <i> <j>              (m lines, 0-based indices)
// Output:
//   SMILES <string>
//   FALLBACK <reason>
#include <iostream>
#include <string>
#include <vector>

#include "smiles.h"

int main() {
    int n = 0;
    int m = 0;
    if (!(std::cin >> n >> m) || n <= 0 || m < 0) {
        std::cerr << "invalid header\n";
        return 2;
    }
    std::vector<reax_smiles::Atom> atoms;
    atoms.reserve(n);
    for (int i = 0; i < n; ++i) {
        std::string element;
        int charge = 0;
        if (!(std::cin >> element)) {
            std::cerr << "missing atom symbol\n";
            return 2;
        }
        // Optional charge is only read when provided on the same token stream.
        // Keep the format simple: element only.
        atoms.push_back({i, element, charge});
    }
    std::vector<reax_smiles::Bond> bonds;
    bonds.reserve(m);
    for (int k = 0; k < m; ++k) {
        int i = 0;
        int j = 0;
        if (!(std::cin >> i >> j)) {
            std::cerr << "missing bond\n";
            return 2;
        }
        bonds.push_back({i, j});
    }

    // Mirror Molecule::update_topology: over-valent N/O get a formal charge
    // so the emitted SMILES is parseable even though the test input has no
    // charge field.
    std::vector<int> degree(n, 0);
    for (const auto& bond : bonds) {
        degree[bond.i]++;
        degree[bond.j]++;
    }
    for (int i = 0; i < n; ++i) {
        if (atoms[i].element == "N" && degree[i] > 3) atoms[i].charge = degree[i] - 3;
        else if (atoms[i].element == "O" && degree[i] > 2) atoms[i].charge = degree[i] - 2;
    }

    reax_smiles::Result result = reax_smiles::canonical_smiles(atoms, bonds);
    if (result.fallback || result.smiles.empty()) {
        std::cout << "FALLBACK " << result.error << "\n";
        return 0;
    }
    std::cout << "SMILES " << result.smiles << "\n";
    return 0;
}
