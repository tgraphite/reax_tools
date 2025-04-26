#include "defines.h"

#include <map>
#include <string>
#include <vector>

int atom_reserve = 10000;
int bond_reserve = 30000;
int mol_reserve = 5000;
int angle_reserve = 60000;
int dihedral_reserve = 100000;

std::vector<std::string> all_elements = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl",
    "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
    "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
    "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
    "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "X"};

std::vector<std::string> sorted_elements = {
    "C",  "H",  "O",  "N",  "S",  "F",  "P",  "He", "Li", "Be", "B",  "Ne", "Na", "Mg", "Al", "Si", "Cl",
    "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
    "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
    "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
    "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "X"};

std::map<std::string, float> default_atomic_radius = {
    {"H", 1.05f},  {"He", 1.4f},  {"Li", 1.82f}, {"Be", 2.f},   {"B", 2.f},    {"C", 1.7f},   {"N", 1.55f},
    {"O", 1.7f},   {"F", 1.35f},  {"Ne", 1.54f}, {"Na", 1.36f}, {"Mg", 1.18f}, {"Al", 2.f},   {"Si", 2.1f},
    {"P", 1.8f},   {"S", 1.8f},   {"Cl", 2.27f}, {"Ar", 1.88f}, {"K", 1.76f},  {"Ca", 1.37f}, {"Sc", 2.f},
    {"Ti", 2.f},   {"V", 2.f},    {"Cr", 2.f},   {"Mn", 2.f},   {"Fe", 2.f},   {"Co", 2.f},   {"Ni", 1.63f},
    {"Cu", 1.4f},  {"Zn", 1.39f}, {"Ga", 1.07f}, {"Ge", 2.f},   {"As", 1.85f}, {"Se", 1.9f},  {"Br", 1.85f},
    {"Kr", 2.02f}, {"Rb", 2.f},   {"Sr", 2.f},   {"Y", 2.f},    {"Zr", 2.f},   {"Nb", 2.f},   {"Mo", 2.f},
    {"Tc", 2.f},   {"Ru", 2.f},   {"Rh", 2.f},   {"Pd", 1.63f}, {"Ag", 1.72f}, {"Cd", 1.58f}, {"In", 1.93f},
    {"Sn", 2.17f}, {"Sb", 2.f},   {"Te", 2.06f}, {"I", 1.98f},  {"Xe", 2.16f}, {"Cs", 2.1f},  {"Ba", 2.f},
    {"La", 2.f},   {"Ce", 2.f},   {"Pr", 2.f},   {"Nd", 2.f},   {"Pm", 2.f},   {"Sm", 2.f},   {"Eu", 2.f},
    {"Gd", 2.f},   {"Tb", 2.f},   {"Dy", 2.f},   {"Ho", 2.f},   {"Er", 2.f},   {"Tm", 2.f},   {"Yb", 2.f},
    {"Lu", 2.f},   {"Hf", 2.f},   {"Ta", 2.f},   {"W", 2.f},    {"Re", 2.f},   {"Os", 2.f},   {"Ir", 2.f},
    {"Pt", 1.72f}, {"Au", 1.66f}, {"Hg", 1.55f}, {"Tl", 1.96f}, {"Pb", 2.02f}, {"Bi", 2.f},   {"Po", 2.f},
    {"At", 2.f},   {"Rn", 2.f},   {"Fr", 2.f},   {"Ra", 2.f},   {"X", 0.0f}};

std::map<std::string, int> default_atom_valence_electrons = {
    {"H", 1},  {"He", 2}, {"Li", 1}, {"Be", 2}, {"B", 3},  {"C", 4},  {"N", 3},  {"O", 2},   {"F", 1},  {"Ne", 8},
    {"Na", 1}, {"Mg", 2}, {"Al", 3}, {"Si", 4}, {"P", 5},  {"S", 6},  {"Cl", 7}, {"Ar", 8},  {"K", 1},  {"Ca", 2},
    {"Sc", 3}, {"Ti", 4}, {"V", 5},  {"Cr", 6}, {"Mn", 7}, {"Fe", 5}, {"Co", 3}, {"Ni", 2},  {"Cu", 2}, {"Zn", 2},
    {"Ga", 8}, {"Ge", 8}, {"As", 5}, {"Se", 8}, {"Br", 7}, {"Kr", 8}, {"Rb", 1}, {"Sr", 2},  {"Y", 3},  {"Zr", 4},
    {"Nb", 5}, {"Mo", 6}, {"Tc", 8}, {"Ru", 8}, {"Rh", 8}, {"Pd", 8}, {"Ag", 1}, {"Cd", 12}, {"In", 3}, {"Sn", 4},
    {"Sb", 8}, {"Te", 8}, {"I", 7},  {"Xe", 8}, {"Cs", 1}, {"Ba", 2}, {"La", 8}, {"Ce", 04}, {"Pr", 8}, {"Nd", 8},
    {"Pm", 8}, {"Sm", 8}, {"Eu", 8}, {"Gd", 8}, {"Tb", 8}, {"Dy", 8}, {"Ho", 8}, {"Er", 8},  {"Tm", 8}, {"Yb", 8},
    {"Lu", 8}, {"Hf", 8}, {"Ta", 8}, {"W", 6},  {"Re", 8}, {"Os", 8}, {"Ir", 8}, {"Pt", 8},  {"Au", 5}, {"Hg", 8},
    {"Tl", 8}, {"Pb", 3}, {"Bi", 5}, {"Po", 8}, {"At", 8}, {"Rn", 8}, {"Fr", 1}, {"Ra", 2},  {"X", 0}};
