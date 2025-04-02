#pragma once

#include <array>
#include <map>
#include <string>
#include <vector>

extern int atom_reserve;
extern int bond_reserve;
extern int mol_reserve;
extern int angle_reserve;
extern int dihedral_reserve;

extern std::vector<std::string> all_elements;
extern std::vector<std::string> sorted_elements;
extern std::map<std::string, float> default_atomic_radius;