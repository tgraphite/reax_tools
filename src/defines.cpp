#include "defines.h"

#include <map>
#include <string>
#include <vector>

int atom_reserve = 10000;
int bond_reserve = 30000;
int mol_reserve = 5000;
int angle_reserve = 60000;
int dihedral_reserve = 100000;
int MAX_RING_SIZE = 8;

std::vector<std::string> all_elements = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl",
    "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
    "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
    "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
    "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "X"};

std::map<std::string, int> element_to_index = {
    {"H", 1},    {"He", 2},   {"Li", 3},   {"Be", 4},   {"B", 5},    {"C", 6},    {"N", 7},    {"O", 8},    {"F", 9},
    {"Ne", 10},  {"Na", 11},  {"Mg", 12},  {"Al", 13},  {"Si", 14},  {"P", 15},   {"S", 16},   {"Cl", 17},  {"Ar", 18},
    {"K", 19},   {"Ca", 20},  {"Sc", 21},  {"Ti", 22},  {"V", 23},   {"Cr", 24},  {"Mn", 25},  {"Fe", 26},  {"Co", 27},
    {"Ni", 28},  {"Cu", 29},  {"Zn", 30},  {"Ga", 31},  {"Ge", 32},  {"As", 33},  {"Se", 34},  {"Br", 35},  {"Kr", 36},
    {"Rb", 37},  {"Sr", 38},  {"Y", 39},   {"Zr", 40},  {"Nb", 41},  {"Mo", 42},  {"Tc", 43},  {"Ru", 44},  {"Rh", 45},
    {"Pd", 46},  {"Ag", 47},  {"Cd", 48},  {"In", 49},  {"Sn", 50},  {"Sb", 51},  {"Te", 52},  {"I", 53},   {"Xe", 54},
    {"Cs", 55},  {"Ba", 56},  {"La", 57},  {"Ce", 58},  {"Pr", 59},  {"Nd", 60},  {"Pm", 61},  {"Sm", 62},  {"Eu", 63},
    {"Gd", 64},  {"Tb", 65},  {"Dy", 66},  {"Ho", 67},  {"Er", 68},  {"Tm", 69},  {"Yb", 70},  {"Lu", 71},  {"Hf", 72},
    {"Ta", 73},  {"W", 74},   {"Re", 75},  {"Os", 76},  {"Ir", 77},  {"Pt", 78},  {"Au", 79},  {"Hg", 80},  {"Tl", 81},
    {"Pb", 82},  {"Bi", 83},  {"Po", 84},  {"At", 85},  {"Rn", 86},  {"Fr", 87},  {"Ra", 88},  {"Ac", 89},  {"Th", 90},
    {"Pa", 91},  {"U", 92},   {"Np", 93},  {"Pu", 94},  {"Am", 95},  {"Cm", 96},  {"Bk", 97},  {"Cf", 98},  {"Es", 99},
    {"Fm", 100}, {"Md", 101}, {"No", 102}, {"Lr", 103}, {"Rf", 104}, {"Db", 105}, {"Sg", 106}, {"Bh", 107}, {"Hs", 108},
    {"Mt", 109}, {"Ds", 110}, {"Rg", 111}, {"Cn", 112}, {"Nh", 113}, {"Fl", 114}, {"Mc", 115}, {"Lv", 116}, {"Ts", 117},
    {"Og", 118}, {"X", 0}};

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
    {"Cu", 1.9f},  {"Zn", 1.39f}, {"Ga", 1.07f}, {"Ge", 2.f},   {"As", 1.85f}, {"Se", 1.9f},  {"Br", 1.85f},
    {"Kr", 2.02f}, {"Rb", 2.f},   {"Sr", 2.f},   {"Y", 2.f},    {"Zr", 2.f},   {"Nb", 2.f},   {"Mo", 2.f},
    {"Tc", 2.f},   {"Ru", 2.f},   {"Rh", 2.f},   {"Pd", 1.63f}, {"Ag", 1.72f}, {"Cd", 1.58f}, {"In", 1.93f},
    {"Sn", 2.17f}, {"Sb", 2.f},   {"Te", 2.06f}, {"I", 1.98f},  {"Xe", 2.16f}, {"Cs", 2.1f},  {"Ba", 2.f},
    {"La", 2.f},   {"Ce", 2.f},   {"Pr", 2.f},   {"Nd", 2.f},   {"Pm", 2.f},   {"Sm", 2.f},   {"Eu", 2.f},
    {"Gd", 2.f},   {"Tb", 2.f},   {"Dy", 2.f},   {"Ho", 2.f},   {"Er", 2.f},   {"Tm", 2.f},   {"Yb", 2.f},
    {"Lu", 2.f},   {"Hf", 2.f},   {"Ta", 2.f},   {"W", 2.f},    {"Re", 2.f},   {"Os", 2.f},   {"Ir", 2.f},
    {"Pt", 1.72f}, {"Au", 1.66f}, {"Hg", 1.55f}, {"Tl", 1.96f}, {"Pb", 2.02f}, {"Bi", 2.f},   {"Po", 2.f},
    {"At", 2.f},   {"Rn", 2.f},   {"Fr", 2.f},   {"Ra", 2.f},   {"X", 0.0f}};

std::map<std::string, int> max_valences = {
    {"H", 1},  {"He", 0}, {"Li", 1}, {"Be", 2}, {"B", 3},  {"C", 4},  {"N", 5},  {"O", 2},  {"F", 1},  {"Ne", 0},
    {"Na", 1}, {"Mg", 2}, {"Al", 3}, {"Si", 4}, {"P", 5},  {"S", 6},  {"Cl", 7}, {"Ar", 0}, {"K", 1},  {"Ca", 2},
    {"Sc", 3}, {"Ti", 4}, {"V", 5},  {"Cr", 6}, {"Mn", 7}, {"Fe", 6}, {"Co", 3}, {"Ni", 2}, {"Cu", 2}, {"Zn", 2},
    {"Ga", 3}, {"Ge", 4}, {"As", 5}, {"Se", 0}, {"Br", 7}, {"Kr", 0}, {"Rb", 1}, {"Sr", 2}, {"Y", 3},  {"Zr", 4},
    {"Nb", 5}, {"Mo", 6}, {"Tc", 8}, {"Ru", 8}, {"Rh", 8}, {"Pd", 8}, {"Ag", 8}, {"Cd", 8}, {"In", 3}, {"Sn", 4},
    {"Sb", 8}, {"Te", 8}, {"I", 7},  {"Xe", 8}, {"Cs", 1}, {"Ba", 2}, {"La", 0}, {"Ce", 6}, {"Pr", 0}, {"Nd", 0},
    {"Pm", 0}, {"Sm", 0}, {"Eu", 0}, {"Gd", 0}, {"Tb", 0}, {"Dy", 0}, {"Ho", 0}, {"Er", 0}, {"Tm", 0}, {"Yb", 0},
    {"Lu", 0}, {"Hf", 0}, {"Ta", 0}, {"W", 0},  {"Re", 0}, {"Os", 0}, {"Ir", 0}, {"Pt", 0}, {"Au", 0}, {"Hg", 0},
    {"Tl", 0}, {"Pb", 0}, {"Bi", 0}, {"Po", 0}, {"At", 0}, {"Rn", 0}, {"Fr", 0}, {"Ra", 0}, {"X", 0}};

std::vector<std::string> default_order = {
    "He", "Li", "Be", "B",  "C",  "H",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl",
    "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
    "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
    "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
    "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "X"};