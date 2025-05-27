#pragma once

#include <array>
#include <map>
#include <string>
#include <vector>

extern int atom_reserve;
extern int bond_reserve;
extern int mol_reserve;

extern std::map<std::string, int> element_to_index;
extern std::vector<std::string> all_elements;
extern std::vector<std::string> sorted_elements;
extern std::map<std::string, float> default_atomic_radius;
extern std::map<std::string, int> max_valences;
extern std::vector<std::string> default_order;

extern int MAX_RING_SIZE;
