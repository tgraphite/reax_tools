#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <fstream>
#include <iostream>

#include "atom.h"
#include "defines.h"
#include "fmt/format.h"
#include "molecule.h"

// TODO: get bond order by graph algorithm
// maybe via 8-electron rule
void draw_molecule(const Molecule& molecule, const std::string& output_path) {
    // Create a new RDKit molecule
    RDKit::RWMol mol;

    // Add atoms
    std::map<int, int> atom_idx_map;  // Maps our atom id to RDKit atom index
    for (const auto& atom : molecule.mol_atoms) {
        int atomic_num = element_to_index[atom->type_name];

        RDKit::Atom rdkit_atom(atomic_num);
        int idx = mol.addAtom(&rdkit_atom);
        atom_idx_map[atom->id] = idx;
    }

    // Add bonds
    for (const auto& bond : molecule.mol_bonds) {
        int begin_idx = atom_idx_map[bond->atom_i->id];
        int end_idx = atom_idx_map[bond->atom_j->id];
        mol.addBond(begin_idx, end_idx, RDKit::Bond::BondType::UNSPECIFIED);
    }

    // Skip sanitization and directly set necessary properties for drawing
    for (auto atom : mol.atoms()) {
        atom->setNoImplicit(true);  // Prevent implicit H addition
        atom->setNumExplicitHs(0);  // Set explicit H count to 0
    }

    // Generate 2D coordinates
    RDKit::DGeomHelpers::EmbedMolecule(mol);

    // Create drawer with custom options
    RDKit::MolDraw2DSVG drawer(400, 400);
    RDKit::MolDraw2DUtils::prepareMolForDrawing(mol);

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
    drawer.drawOptions().prepareMolsBeforeDrawing = false;

    // Draw the molecule
    drawer.drawMolecule(mol);
    drawer.finishDrawing();

    // Save to file
    std::ofstream out(output_path);
    out << drawer.getDrawingText();
    out.close();
}
