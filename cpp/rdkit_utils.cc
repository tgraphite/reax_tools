#include "rdkit_utils.h"

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <filesystem>
#include <string>

#include "atom.h"
#include "fmt/format.h"
#include "molecule.h"

// Main SMILES generator, RDKit-style
std::string rdkit_smiles(const Molecule& mol) {
    RDKit::RWMol rdkit_mol;
    std::map<int, int> atom_idx_map;  // Maps our atom id to RDKit atom index
    for (const auto& atom : mol.mol_atoms) {
        int atomic_num = element_to_index[atom->type_name];

        RDKit::Atom rdkit_atom(atomic_num);
        int idx = rdkit_mol.addAtom(&rdkit_atom);
        atom_idx_map[atom->id] = idx;
    }

    // Add bonds
    for (const auto& bond : mol.mol_bonds) {
        int begin_idx = atom_idx_map[bond->atom_i->id];
        int end_idx = atom_idx_map[bond->atom_j->id];

        // If bond already exists (rarely happens), skip
        if (rdkit_mol.getBondBetweenAtoms(begin_idx, end_idx)) {
            continue;
        }

        if (bond->order == 1) {
            rdkit_mol.addBond(begin_idx, end_idx, RDKit::Bond::BondType::SINGLE);
        } else if (bond->order == 2) {
            rdkit_mol.addBond(begin_idx, end_idx, RDKit::Bond::BondType::DOUBLE);
        } else if (bond->order == 3) {
            rdkit_mol.addBond(begin_idx, end_idx, RDKit::Bond::BondType::TRIPLE);
        } else {
            rdkit_mol.addBond(begin_idx, end_idx, RDKit::Bond::BondType::UNSPECIFIED);
        }

        // mol.addBond(begin_idx, end_idx, RDKit::Bond::BondType::UNSPECIFIED);
    }

    // Skip sanitization and directly set necessary properties for drawing
    for (auto atom : rdkit_mol.atoms()) {
        atom->setNoImplicit(true);  // Prevent implicit H addition
        atom->setNumExplicitHs(0);  // Set explicit H count to 0
    }

    std::string smiles = RDKit::MolToSmiles(rdkit_mol);
    return smiles;
}

// Main SMILES generator, RDKit-style
void rdkit_draw_molecule(const Molecule& mol, std::string output_dir) {
    // Too big, don't draw, just use formula in downstream visualization tool.
    if (mol.mol_atoms.size() > 60) return;

    RDKit::RWMol rdkit_mol;
    std::map<int, int> atom_idx_map;  // Maps our atom id to RDKit atom index
    for (const auto& atom : mol.mol_atoms) {
        int atomic_num = element_to_index[atom->type_name];

        RDKit::Atom rdkit_atom(atomic_num);
        int idx = rdkit_mol.addAtom(&rdkit_atom);
        atom_idx_map[atom->id] = idx;
    }

    // Add bonds
    for (const auto& bond : mol.mol_bonds) {
        int begin_idx = atom_idx_map[bond->atom_i->id];
        int end_idx = atom_idx_map[bond->atom_j->id];

        // If bond already exists (rarely happens), skip
        if (rdkit_mol.getBondBetweenAtoms(begin_idx, end_idx)) {
            continue;
        }

        if (bond->order == 1) {
            rdkit_mol.addBond(begin_idx, end_idx, RDKit::Bond::BondType::SINGLE);
        } else if (bond->order == 2) {
            rdkit_mol.addBond(begin_idx, end_idx, RDKit::Bond::BondType::DOUBLE);
        } else if (bond->order == 3) {
            rdkit_mol.addBond(begin_idx, end_idx, RDKit::Bond::BondType::TRIPLE);
        } else {
            rdkit_mol.addBond(begin_idx, end_idx, RDKit::Bond::BondType::UNSPECIFIED);
        }

        // mol.addBond(begin_idx, end_idx, RDKit::Bond::BondType::UNSPECIFIED);
    }

    // Skip sanitization and directly set necessary properties for drawing
    for (auto atom : rdkit_mol.atoms()) {
        atom->setNoImplicit(true);  // Prevent implicit H addition
        atom->setNumExplicitHs(0);  // Set explicit H count to 0
    }

    RDKit::DGeomHelpers::EmbedMolecule(rdkit_mol);

    // Create drawer with custom options
    RDKit::MolDraw2DSVG drawer(250, 250);
    RDKit::MolDraw2DUtils::prepareMolForDrawing(rdkit_mol);

    // Configure drawer options
    drawer.drawOptions().addAtomIndices = false;
    drawer.drawOptions().addBondIndices = false;
    drawer.drawOptions().includeAtomTags = false;
    drawer.drawOptions().includeRadicals = true;
    drawer.drawOptions().atomLabelDeuteriumTritium = false;
    drawer.drawOptions().explicitMethyl = false;
    drawer.drawOptions().includeChiralFlagLabel = false;
    drawer.drawOptions().comicMode = false;
    drawer.drawOptions().continuousHighlight = false;
    drawer.drawOptions().circleAtoms = false;
    drawer.drawOptions().dummiesAreAttachments = false;
    drawer.drawOptions().includeMetadata = false;
    drawer.drawOptions().prepareMolsBeforeDrawing = true;
    drawer.drawOptions().centreMoleculesBeforeDrawing = true;
    // Draw the molecule
    drawer.drawMolecule(rdkit_mol);
    drawer.finishDrawing();

    // Create output directory if it doesn't exist
    std::string pictures_dir = output_dir + "molecule_pictures/";
    std::filesystem::create_directories(pictures_dir);

    // Save to file
    std::string output_path = output_dir + fmt::format("molecule_pictures/molecule_{}.svg", mol.formula);
    std::ofstream out(output_path);
    out << drawer.getDrawingText();
    out.close();
}