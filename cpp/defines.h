#pragma once

#include <map>
#include <string>
#include <vector>

extern int atom_reserve;
extern int bond_reserve;
extern int mol_reserve;
extern int angle_reserve;
extern int dihedral_reserve;
extern int MAX_RING_SIZE;

extern std::vector<std::string> all_elements;
extern std::map<std::string, int> element_to_index;
extern std::map<std::string, float> default_atomic_radius;
extern std::map<std::string, int> max_valences;
extern std::vector<std::string> default_order;
extern std::vector<unsigned int> prime_numbers;
extern std::vector<unsigned int> bigger_prime_numbers;
