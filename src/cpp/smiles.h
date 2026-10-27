#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace reax_smiles {

// Minimal explicit-atom graph. Bonds are connectivity-only in this version:
// bond order is deliberately not part of the structural identity.
struct Atom {
    int id = 0;
    std::string element;  // "C", "H", "Cl", ...; unknown -> "*"
    int charge = 0;
};

struct Bond {
    int i = 0;  // 0-based index into the atom vector
    int j = 0;
};

// Components above this size are formula-level identities.
constexpr std::size_t kMaxSmilesAtoms = 200;

struct Options {
    // Molecules above this size are treated as formula-level identities:
    // no SMILES is generated and the caller uses the formula hash. 200 is a
    // deliberate cutoff for kerogen/polymer-scale components, where the
    // probability of a meaningful constitutional isomer pair is negligible.
    std::size_t max_atoms = kMaxSmilesAtoms;
    // Colour refinement rounds before giving up.
    std::size_t max_refine_rounds = 4096;
    // Maximum explored IR leaves before giving up.
    std::size_t max_ir_leaves = 20000;
    // Maximum size of an individualisation cell before giving up.
    std::size_t max_cell = 512;
};

struct Result {
    std::string smiles;     // empty when fallback/error
    bool fallback = false;  // oversized or search budget exceeded
    std::string error;      // diagnostic only
};

// Canonical SMILES with every hydrogen written explicitly as [H].
// The graph must contain all atoms, including hydrogens.
Result canonical_smiles(const std::vector<Atom>& atoms,
                        const std::vector<Bond>& bonds,
                        const Options& options = {});

// Deterministic 64-bit FNV-1a digests.  The domain prefix prevents a
// formula string and a SMILES string from sharing an ID by construction.
std::uint64_t smiles_id(const std::string& smiles);
std::uint64_t formula_id(const std::string& formula);

}  // namespace reax_smiles
