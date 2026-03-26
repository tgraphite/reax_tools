#include "argparser.h"

#include <map>
#include <sstream>
#include <string>
#include <vector>

int MAX_RING_SIZE = 8;
int MIN_RING_SIZE = 5;
int MAX_NEIGH = 10;

int NUM_THREADS = 4;
int MAX_REACTIONS = 60;
int MAX_KEY_MOLECULES = 60;
int DUMP_STEPS = 10;
int NETWORK_FLOW_MAX_REACTANTS = 10;
int NETWORK_FLOW_MAX_PRODUCTS = 10;

std::string INPUT_FILE = "";
std::string OUTPUT_DIR = "reax_tools_output/";
std::vector<std::string> INPUT_ELEMENT_TYPES = {};
float RVDW_FACTOR = 1.2;

bool FLAG_DUMP_STRUCTURE = false;
bool FLAG_MARK_RING_ATOMS = false;

bool FLAG_NO_REDUCE_REACTIONS = false;
bool FLAG_NO_RINGS = false;
bool FLAG_NO_REACTIONS = false;

// Reaction tracking flags (ChemTrayzer-like)
bool FLAG_TRACK_REACTIONS = true;
bool FLAG_NO_TRACK_REACTIONS = false;
int STABLE_TIME_FRAMES = 3;
float TIMESTEP_FS = 0.25f;
int SAMPLING_FREQ = 10;

bool FLAG_MERGE_ELEMENTS = false;
bool FLAG_RESCALE_MERGE_COUNT = false;
std::string MERGE_TARGET = "";
std::vector<int> MERGE_RANGES = {};

std::vector<std::string> ALL_ELEMENTS = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl",
    "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
    "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
    "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
    "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "X"};

std::map<std::string, int> ELEMENT_TO_INDEX = {
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

std::map<std::string, float> ELEMENT_ATOMIC_RADII = {
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

std::map<std::string, int> ELEMENT_MAX_VALENCIES = {
    {"H", 1},  {"He", 0}, {"Li", 1}, {"Be", 2}, {"B", 3},  {"C", 4},  {"N", 4},  {"O", 2},  {"F", 1},  {"Ne", 0},
    {"Na", 1}, {"Mg", 2}, {"Al", 3}, {"Si", 4}, {"P", 6},  {"S", 6},  {"Cl", 7}, {"Ar", 0}, {"K", 1},  {"Ca", 2},
    {"Sc", 3}, {"Ti", 4}, {"V", 5},  {"Cr", 6}, {"Mn", 7}, {"Fe", 6}, {"Co", 3}, {"Ni", 2}, {"Cu", 2}, {"Zn", 2},
    {"Ga", 3}, {"Ge", 4}, {"As", 5}, {"Se", 0}, {"Br", 7}, {"Kr", 0}, {"Rb", 1}, {"Sr", 2}, {"Y", 3},  {"Zr", 4},
    {"Nb", 5}, {"Mo", 6}, {"Tc", 8}, {"Ru", 8}, {"Rh", 8}, {"Pd", 8}, {"Ag", 8}, {"Cd", 8}, {"In", 3}, {"Sn", 4},
    {"Sb", 8}, {"Te", 8}, {"I", 7},  {"Xe", 8}, {"Cs", 1}, {"Ba", 2}, {"La", 0}, {"Ce", 6}, {"Pr", 0}, {"Nd", 0},
    {"Pm", 0}, {"Sm", 0}, {"Eu", 0}, {"Gd", 0}, {"Tb", 0}, {"Dy", 0}, {"Ho", 0}, {"Er", 0}, {"Tm", 0}, {"Yb", 0},
    {"Lu", 0}, {"Hf", 0}, {"Ta", 0}, {"W", 0},  {"Re", 0}, {"Os", 0}, {"Ir", 0}, {"Pt", 0}, {"Au", 0}, {"Hg", 0},
    {"Tl", 0}, {"Pb", 0}, {"Bi", 0}, {"Po", 0}, {"At", 0}, {"Rn", 0}, {"Fr", 0}, {"Ra", 0}, {"X", 0}};

std::vector<std::string> ELEMENT_DISPLAY_ORDER = {
    "He", "Li", "Be", "B",  "C",  "H",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl",
    "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
    "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
    "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
    "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "X"};

std::vector<unsigned int> PRIME_NUMBERS = {
    3,   5,   7,   11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,  61,  67,  71,  73,
    79,  83,  89,  97,  101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
    181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283,
    293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419,
    421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
    557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659};

std::vector<unsigned int> BIGGER_PRIME_NUMBERS = {
    661,  673,  677,  683,  691,  701,  709,  719,  727,  733,  739,  743,  751,  757,  761,  769,  773,
    787,  797,  809,  811,  821,  823,  827,  829,  839,  853,  857,  859,  863,  877,  881,  883,  887,
    907,  911,  919,  929,  937,  941,  947,  953,  967,  971,  977,  983,  991,  997,  1009, 1013, 1019,
    1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123,
    1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259,
    1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399,
    1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499};

std::map<std::string, float> ELEMENT_MASS = {
    {"H", 1.0080},      {"He", 4.00260},   {"Li", 7.0},        {"Be", 9.012183},   {"B", 10.81},
    {"C", 12.011},      {"N", 14.007},     {"O", 15.999},      {"F", 18.99840316}, {"Ne", 20.180},
    {"Na", 22.9897693}, {"Mg", 24.305},    {"Al", 26.981538},  {"Si", 28.085},     {"P", 30.97376200},
    {"S", 32.07},       {"Cl", 35.45},     {"Ar", 39.9},       {"K", 39.0983},     {"Ca", 40.08},
    {"Sc", 44.95591},   {"Ti", 47.867},    {"V", 50.9415},     {"Cr", 51.996},     {"Mn", 54.93804},
    {"Fe", 55.84},      {"Co", 58.93319},  {"Ni", 58.693},     {"Cu", 63.55},      {"Zn", 65.4},
    {"Ga", 69.723},     {"Ge", 72.63},     {"As", 74.92159},   {"Se", 78.97},      {"Br", 79.90},
    {"Kr", 83.80},      {"Rb", 85.468},    {"Sr", 87.62},      {"Y", 88.90584},    {"Zr", 91.22},
    {"Nb", 92.90637},   {"Mo", 95.95},     {"Tc", 96.90636},   {"Ru", 101.1},      {"Rh", 102.9055},
    {"Pd", 106.42},     {"Ag", 107.868},   {"Cd", 112.41},     {"In", 114.818},    {"Sn", 118.71},
    {"Sb", 121.760},    {"Te", 127.6},     {"I", 126.9045},    {"Xe", 131.29},     {"Cs", 132.9054520},
    {"Ba", 137.33},     {"La", 138.9055},  {"Ce", 140.116},    {"Pr", 140.90766},  {"Nd", 144.24},
    {"Pm", 144.91276},  {"Sm", 150.4},     {"Eu", 151.964},    {"Gd", 157.25},     {"Tb", 158.92535},
    {"Dy", 162.500},    {"Ho", 164.93033}, {"Er", 167.26},     {"Tm", 168.93422},  {"Yb", 173.05},
    {"Lu", 174.9667},   {"Hf", 178.49},    {"Ta", 180.9479},   {"W", 183.84},      {"Re", 186.207},
    {"Os", 190.2},      {"Ir", 192.22},    {"Pt", 195.08},     {"Au", 196.96657},  {"Hg", 200.59},
    {"Tl", 204.383},    {"Pb", 207},       {"Bi", 208.98040},  {"Po", 208.98243},  {"At", 209.98715},
    {"Rn", 222.01758},  {"Fr", 223.01973}, {"Ra", 226.02541},  {"Ac", 227.02775},  {"Th", 232.038},
    {"Pa", 231.03588},  {"U", 238.0289},   {"Np", 237.048172}, {"Pu", 244.06420},  {"Am", 243.061380},
    {"Cm", 247.07035},  {"Bk", 247.07031}, {"Cf", 251.07959},  {"Es", 252.0830},   {"Fm", 257.09511},
    {"Md", 258.09843},  {"No", 259.10100}, {"Lr", 266.120},    {"Rf", 267.122},    {"Db", 268.126},
    {"Sg", 269.128},    {"Bh", 270.133},   {"Hs", 269.1336},   {"Mt", 277.154},    {"Ds", 282.166},
    {"Rg", 282.169},    {"Cn", 286.179},   {"Nh", 286.182},    {"Fl", 290.192},    {"Mc", 290.196},
    {"Lv", 293.205},    {"Ts", 294.211},   {"Og", 295.216}};

bool operation_set_file(std::vector<std::string> input_string) {
    if (input_string.size() != 1) return false;
    INPUT_FILE = input_string[0];
    return true;
}
bool operation_set_output_dir(std::vector<std::string> input_string) {
    if (input_string.size() != 1) return false;
    OUTPUT_DIR = input_string[0];
    return true;
}
bool operation_set_rvdw_factor(std::vector<std::string> input_string) {
    if (input_string.size() != 1) return false;
    RVDW_FACTOR = std::stof(input_string[0]);
    return true;
}
bool operation_set_type_radii(std::vector<std::string> input_string) {
    for (const auto& tr : input_string) {
        auto pos = tr.find(':');
        if (pos != std::string::npos) {
            std::string elem = tr.substr(0, pos);
            float radius = std::stof(tr.substr(pos + 1));
            ELEMENT_ATOMIC_RADII[elem] = radius;
        } else {
            return false;
        }
    }
    return true;
}
bool operation_set_types(std::vector<std::string> input_string) {
    INPUT_ELEMENT_TYPES.clear();
    if (input_string.size() == 1) {
        std::stringstream ss(input_string[0]);
        std::string item;
        while (std::getline(ss, item, ',')) {
            INPUT_ELEMENT_TYPES.push_back(item);
        }
    } else {
        for (const auto& s : input_string) {
            INPUT_ELEMENT_TYPES.push_back(s);
        }
    }
    return true;
}
bool operation_set_type_valence(std::vector<std::string> input_string) {
    for (const auto& tv : input_string) {
        auto pos = tv.find(':');
        if (pos != std::string::npos) {
            std::string elem = tv.substr(0, pos);
            int valence = std::stoi(tv.substr(pos + 1));
            ELEMENT_MAX_VALENCIES[elem] = valence;
        } else {
            return false;
        }
    }
    return true;
}
bool operation_set_num_threads(std::vector<std::string> input_string) {
    if (input_string.size() != 1) return false;
    NUM_THREADS = std::stoi(input_string[0]);
    return true;
}
bool operation_set_flag_dump_structure(std::vector<std::string> input_string) {
    if (input_string.size() != 1) return false;
    FLAG_DUMP_STRUCTURE = true;
    DUMP_STEPS = std::stoi(input_string[0]);
    return true;
}
bool operation_set_flag_mark_ring_atoms(std::vector<std::string> input_string) {
    FLAG_MARK_RING_ATOMS = true;
    return true;
}
bool operation_set_flag_no_rings(std::vector<std::string> input_string) {
    FLAG_NO_RINGS = true;
    return true;
}
bool operation_set_flag_no_reactions(std::vector<std::string> input_string) {
    FLAG_NO_REACTIONS = true;
    return true;
}
bool operation_set_flag_no_reduce_reactions(std::vector<std::string> input_string) {
    FLAG_NO_REDUCE_REACTIONS = true;
    return true;
}
bool operation_set_max_reactions(std::vector<std::string> input_string) {
    if (input_string.size() != 1) return false;
    MAX_REACTIONS = std::stoi(input_string[0]);
    return true;
}
bool operation_set_merge_target(std::vector<std::string> input_string) {
    if (input_string.size() != 1) return false;
    MERGE_TARGET = input_string[0];
    FLAG_MERGE_ELEMENTS = true;
    return true;
}
bool operation_set_merge_ranges(std::vector<std::string> input_string) {
    MERGE_RANGES.clear();
    if (input_string.size() == 1) {
        std::stringstream ss(input_string[0]);
        std::string item;
        while (std::getline(ss, item, ',')) {
            MERGE_RANGES.push_back(std::stoi(item));
        }
    } else {
        for (const auto& s : input_string) {
            MERGE_RANGES.push_back(std::stoi(s));
        }
    }
    return true;
}
bool operation_set_rescale_merge_counts(std::vector<std::string> input_string) {
    FLAG_RESCALE_MERGE_COUNT = true;
    return true;
}
bool operation_set_flag_no_track_reactions(std::vector<std::string> input_string) {
    FLAG_TRACK_REACTIONS = false;
    FLAG_NO_TRACK_REACTIONS = true;
    return true;
}
bool operation_set_stable_time(std::vector<std::string> input_string) {
    if (input_string.size() != 1) return false;
    STABLE_TIME_FRAMES = std::stoi(input_string[0]);
    return true;
}
bool operation_set_timestep(std::vector<std::string> input_string) {
    if (input_string.size() != 1) return false;
    TIMESTEP_FS = std::stof(input_string[0]);
    return true;
}
bool operation_set_sampling_freq(std::vector<std::string> input_string) {
    if (input_string.size() != 1) return false;
    SAMPLING_FREQ = std::stoi(input_string[0]);
    return true;
}
bool operation_set_element_order(std::vector<std::string> input_string) {
    ELEMENT_DISPLAY_ORDER.clear();
    if (input_string.size() == 1) {
        std::stringstream ss(input_string[0]);
        std::string item;
        while (std::getline(ss, item, ',')) {
            ELEMENT_DISPLAY_ORDER.push_back(item);
        }
    } else {
        for (const auto& s : input_string) {
            ELEMENT_DISPLAY_ORDER.push_back(s);
        }
    }
    return true;
}

ArgParser init_argparser() {
    ArgParser parser("reax_tools", "ReaxTools: High performance reactive MD post-processer");

    parser.add_argument("--traj", "-f", "input trajectory file (.xyz/.lammpstrj)", "key input", false,
                        operation_set_file);
    parser.add_argument("--output", "-o", "output directory", "key input", false, operation_set_output_dir);
    parser.add_argument("--rescale-vdw", "-r", "rescale vdw radii by factor (default=1.2)", "key input", false,
                        operation_set_rvdw_factor);
    parser.add_argument("--types", "-t", "element types splitted in commas", "key input", false, operation_set_types);
    parser.add_argument("--type-radius", "-tr", "set atomic radius for element, e.g. -tr N:1.5", "key input", false,
                        operation_set_type_radii);
    parser.add_argument("--type-valence", "-tv", "set atomic valence for element, e.g. -tv N:4", "key input", false,
                        operation_set_type_valence);
    parser.add_argument("--threads", "-nt", "number of threads", "Performance", false, operation_set_num_threads);
    parser.add_argument("--dump", "", "dump lammps data file (.data) for every N frames", "Misc", false,
                        operation_set_flag_dump_structure);
    parser.add_argument("--mark-ring-atoms", "", "use mol id entry to mark ring atom when dump data", "Misc", true,
                        operation_set_flag_mark_ring_atoms);
    parser.add_argument("--no-rings", "", "disable ring computation", "Misc", true, operation_set_flag_no_rings);
    parser.add_argument("--no-reactions", "", "disable reaction analysis", "Network output options", true,
                        operation_set_flag_no_reactions);
    parser.add_argument("--no-reduce-reactions", "-norr", "disable reduce reverse reactions", "Network output options",
                        true, operation_set_flag_no_reduce_reactions);
    parser.add_argument("--max-reactions", "", "set max reactions in default network", "Network output options", false,
                        operation_set_max_reactions);
    parser.add_argument("--merge-element", "-me", "merge species groups by an element type", "Species output options",
                        false, operation_set_merge_target);
    parser.add_argument("--merge-ranges", "-mr", "merge group range, split by commas", "Species output options", false,
                        operation_set_merge_ranges);
    parser.add_argument("--rescale-count", "-rc",
                        "rescale group weight by actual numbers of atoms (instead of molecules)",
                        "Species output options", true, operation_set_rescale_merge_counts);
    parser.add_argument("--element-order", "-order", "set element order of outputting formulas",
                        "Species output options", false, operation_set_element_order);
    parser.add_argument("--no-track-reactions", "", "disable ChemTrayzer-like reaction tracking", "Reaction tracking",
                        true, operation_set_flag_no_track_reactions);
    parser.add_argument("--stable-time", "", "minimum frames for molecule to be stable (default=3)", "Reaction tracking",
                        false, operation_set_stable_time);
    parser.add_argument("--timestep", "", "MD timestep in fs (default=0.25)", "Reaction tracking", false,
                        operation_set_timestep);
    parser.add_argument("--sampling-freq", "", "frame sampling frequency (default=10)", "Reaction tracking", false,
                        operation_set_sampling_freq);

    return parser;
}