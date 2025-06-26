#pragma once

#include <string>

#include "system.h"

std::string rdkit_smiles(const Molecule &mol);
void rdkit_draw_molecule(const Molecule &mol, std::string output_dir);