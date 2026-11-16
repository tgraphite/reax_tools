#include "system.h"

#include <assert.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <mutex>
#include <queue>
#include <sstream>

#include "argparser.h"
#include "cell_list.h"
#include "fmt/format.h"
#include "reax_counter.h"
#include "smiles.h"
#include "universe.h"
#include "vec_algorithms.h"

std::mutex species_counter_mutex;

namespace {

// Type ids must be stable for the whole trajectory.  Otherwise a frame whose
// first atom happens to be O assigns type 1 to O, and a later frame whose first
// atom is H reassigns type 1 to H, so a sequence of --dump files no longer
// agrees on element names.
std::map<std::string, int> global_type_stoi;
std::map<int, std::string> global_type_itos;
bool global_type_registry_initialized = false;

struct PropertyColumn {
    std::string name;
    char type = 'R';
    int count = 1;
    int offset = 0;
};

bool is_integer_token(const std::string& s) {
    if (s.empty()) return false;
    for (char ch : s) {
        if (!std::isdigit(static_cast<unsigned char>(ch))) return false;
    }
    return true;
}

bool is_number_token(const std::string& s) {
    if (s.empty()) return false;
    char* end = nullptr;
    std::strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0';
}

// Parse an extended-XYZ comment line into key=value pairs.  Quoted values may
// contain spaces; bracket values may contain commas/spaces.  Unknown keys are
// ignored by the caller.
void parse_key_values(const std::string& line, std::map<std::string, std::string>& out) {
    std::size_t i = 0;
    while (i < line.size()) {
        while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
        if (i >= line.size() || line[i] == '#') break;

        const std::size_t key_start = i;
        while (i < line.size() && line[i] != '=' && !std::isspace(static_cast<unsigned char>(line[i]))) ++i;
        std::string key = line.substr(key_start, i - key_start);
        while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
        if (i >= line.size() || line[i] != '=') {
            while (i < line.size() && !std::isspace(static_cast<unsigned char>(line[i]))) ++i;
            continue;
        }
        ++i;
        while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;

        std::string value;
        if (i < line.size() && line[i] == '"') {
            ++i;
            while (i < line.size()) {
                const char ch = line[i];
                if (ch == '\\' && i + 1 < line.size()) {
                    value += line[i + 1];
                    i += 2;
                    continue;
                }
                if (ch == '"') {
                    ++i;
                    break;
                }
                value += ch;
                ++i;
            }
        }
        else if (i < line.size() && line[i] == '[') {
            int depth = 0;
            while (i < line.size()) {
                const char ch = line[i];
                value += ch;
                if (ch == '[') ++depth;
                if (ch == ']') {
                    --depth;
                    if (depth == 0) {
                        ++i;
                        break;
                    }
                }
                ++i;
            }
        }
        else {
            while (i < line.size() && !std::isspace(static_cast<unsigned char>(line[i]))) {
                value += line[i];
                ++i;
            }
        }
        if (!key.empty()) out[key] = value;
    }
}

std::vector<std::string> list_tokens(std::string value) {
    for (char& ch : value) {
        if (ch == '[' || ch == ']' || ch == ',') ch = ' ';
    }
    return split_by_space(value);
}

bool parse_double_list(const std::string& value, std::vector<double>& out) {
    out.clear();
    for (const auto& token : list_tokens(value)) {
        try {
            out.push_back(std::stod(token));
        }
        catch (...) {
            return false;
        }
    }
    return !out.empty();
}

bool parse_bool_list(const std::string& value, std::array<bool, 3>& out) {
    const auto tokens = list_tokens(value);
    if (tokens.size() != 3) return false;
    for (int i = 0; i < 3; ++i) {
        std::string t = tokens[i];
        for (char& ch : t) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        if (t == "T" || t == "TRUE" || t == "1") out[i] = true;
        else if (t == "F" || t == "FALSE" || t == "0") out[i] = false;
        else return false;
    }
    return true;
}

std::vector<PropertyColumn> parse_properties(const std::string& value) {
    std::vector<PropertyColumn> columns;
    const auto parts = split(value, ":");
    int offset = 0;
    for (std::size_t i = 0; i + 2 < parts.size(); i += 3) {
        PropertyColumn column;
        column.name = parts[i];
        column.type = parts[i + 1].empty() ? 'R' : parts[i + 1][0];
        try {
            column.count = std::max(1, std::stoi(parts[i + 2]));
        }
        catch (...) {
            column.count = 1;
        }
        column.offset = offset;
        offset += column.count;
        columns.push_back(column);
    }
    return columns;
}

int find_column(const std::vector<PropertyColumn>& columns, std::initializer_list<const char*> names) {
    for (const char* name : names) {
        for (const auto& column : columns) {
            if (column.name == name) return column.offset;
        }
    }
    return -1;
}

int find_column_count(const std::vector<PropertyColumn>& columns, std::initializer_list<const char*> names) {
    for (const char* name : names) {
        for (const auto& column : columns) {
            if (column.name == name) return column.count;
        }
    }
    return 0;
}

bool looks_like_atom_line(const std::string& line) {
    if (line.empty() || line[0] == '#' || line.find('=') != std::string::npos) return false;
    const auto tokens = split_by_space(line);
    if (tokens.size() < 4) return false;
    return is_number_token(tokens[1]) && is_number_token(tokens[2]) && is_number_token(tokens[3]);
}

}  // namespace

System::System() { set_types(); }

System::~System() {
    // No need to free atoms ptr and bonds ptr, they will be deleted in Molecule destrcutor

    for (const auto& instance : molecules) {
        if (instance != nullptr) delete instance;
    }
    molecules.clear();
}

void System::load_xyz(std::ifstream& file) {
    std::string line;
    total_atoms = 0;

    // Find the atom-count line.  A frame is read per call, so the stream is
    // left exactly after this frame's atoms.
    while (std::getline(file, line)) {
        const auto tokens = split_by_space(line);
        if (tokens.empty() || tokens[0] == "#") continue;
        if (tokens.size() == 1 && is_integer_token(tokens[0])) {
            total_atoms = std::stoi(tokens[0]);
            break;
        }
    }
    if (total_atoms <= 0) return;

    atoms.reserve(total_atoms);
    bonds.reserve(total_atoms * 3);
    molecules.reserve(total_atoms / 2);

    std::string comment;
    if (!std::getline(file, comment)) return;

    // A few legacy XYZ files omit the comment line.  Detect that case and
    // treat the line as the first atom rather than silently dropping it.
    bool pending_atom_line = false;
    std::string first_atom_line;
    if (looks_like_atom_line(comment)) {
        pending_atom_line = true;
        first_atom_line = comment;
        comment.clear();
    }

    std::map<std::string, std::string> header;
    parse_key_values(comment, header);

    std::array<double, 3> origin{{0.0, 0.0, 0.0}};
    auto origin_it = header.find("Origin");
    if (origin_it != header.end()) {
        std::vector<double> values;
        if (parse_double_list(origin_it->second, values) && values.size() >= 3) {
            origin = {values[0], values[1], values[2]};
        }
    }

    Cell frame_cell;
    std::array<bool, 3> pbc{{false, false, false}};
    bool lattice_present = false;

    auto lattice_it = header.find("Lattice");
    if (lattice_it != header.end()) {
        std::vector<double> values;
        if (parse_double_list(lattice_it->second, values)) {
            if (values.size() == 3) {
                frame_cell.set_lattice_diagonal({values[0], values[1], values[2]}, origin, {true, true, true});
                lattice_present = true;
            }
            else if (values.size() == 9) {
                // OVITO/ASE write the 9-value form column-major.  The optional
                // 2D form is row-major, so transpose it to the same storage.
                const bool two_dimensional =
                    std::count(lattice_it->second.begin(), lattice_it->second.end(), '[') >= 2;
                if (two_dimensional) {
                    frame_cell.set_vectors({values[0], values[3], values[6]},
                        {values[1], values[4], values[7]}, {values[2], values[5], values[8]}, origin,
                        {true, true, true});
                }
                else {
                    frame_cell.set_lattice9({values[0], values[1], values[2], values[3], values[4], values[5],
                                                values[6], values[7], values[8]},
                        origin, {true, true, true});
                }
                lattice_present = true;
            }
            else {
                fmt::print(stderr, "Warning: Lattice must contain 3 or 9 values, got {}\n", values.size());
            }
        }
        else {
            fmt::print(stderr, "Warning: could not parse Lattice value\n");
        }
    }

    if (lattice_present) {
        pbc = {true, true, true};
        auto pbc_it = header.find("pbc");
        if (pbc_it != header.end()) {
            std::array<bool, 3> parsed{{true, true, true}};
            if (parse_bool_list(pbc_it->second, parsed)) pbc = parsed;
            else fmt::print(stderr, "Warning: could not parse pbc value\n");
        }
        frame_cell.periodic = pbc;
    }
    else {
        frame_cell.origin = origin;
        frame_cell.periodic = {false, false, false};
        frame_cell.valid = false;
    }

    cell = frame_cell;
    has_boundaries = cell.valid && cell.any_periodic();
    if (lattice_present && !cell.valid) {
        fmt::print(stderr, "Warning: singular Lattice, disabling periodic boundaries\n");
    }

    std::vector<PropertyColumn> columns;
    auto properties_it = header.find("Properties");
    if (properties_it != header.end()) columns = parse_properties(properties_it->second);
    if (columns.empty()) {
        columns = {{"species", 'S', 1, 0}, {"pos", 'R', 3, 1}};
    }

    const int species_column = find_column(columns, {"species", "element", "type"});
    const int pos_column = find_column(columns, {"pos", "positions"});
    const int id_column = find_column(columns, {"id"});
    if (species_column < 0 || pos_column < 0 ||
        find_column_count(columns, {"pos", "positions"}) != 3) {
        fmt::print(stderr, "Error: extended XYZ requires species and pos properties\n");
        return;
    }

    int required_columns = 0;
    for (const auto& column : columns) {
        required_columns = std::max(required_columns, column.offset + column.count);
    }

    auto add_type = [&](const std::string& name) -> int {
        auto found = type_stoi.find(name);
        if (found != type_stoi.end()) return found->second;
        const int type_id = System::register_type(name);
        type_itos[type_id] = name;
        type_stoi[name] = type_id;
        total_types = std::max(total_types, type_id);
        return type_id;
    };

    int atom_count = 0;
    while (atom_count < total_atoms) {
        if (pending_atom_line) {
            line = first_atom_line;
            pending_atom_line = false;
        }
        else if (!std::getline(file, line)) {
            break;
        }

        const auto tokens = split_by_space(line);
        if (tokens.empty() || tokens[0] == "#") continue;
        if (static_cast<int>(tokens.size()) < required_columns) {
            fmt::print(stderr, "Error: XYZ atom line has {} columns, expected at least {}\n", tokens.size(),
                required_columns);
            return;
        }

        const std::string raw_type = tokens[species_column];
        int type_id = 0;
        std::string type_name;
        if (is_integer_token(raw_type)) {
            type_id = std::stoi(raw_type);
            auto found = type_itos.find(type_id);
            if (found == type_itos.end()) {
                fmt::print(stderr,
                    "Error: numeric XYZ type {} has no element mapping; provide --types/-t\n", type_id);
                return;
            }
            type_name = found->second;
        }
        else {
            type_name = raw_type;
            type_id = add_type(type_name);
        }

        const int atom_id = (id_column >= 0) ? std::stoi(tokens[id_column]) : atom_count + 1;
        const double x = std::stod(tokens[pos_column]);
        const double y = std::stod(tokens[pos_column + 1]);
        const double z = std::stod(tokens[pos_column + 2]);

        std::array<double, 3> relative{{x - origin[0], y - origin[1], z - origin[2]}};
        if (cell.valid) relative = cell.wrap_relative(relative);

        Atom* atom = new Atom(atom_id, type_id,
            {static_cast<float>(relative[0]), static_cast<float>(relative[1]), static_cast<float>(relative[2])},
            type_name);
        atoms.push_back(atom);
        id_atom_map[atom_id] = atom;
        ++atom_count;
    }

    if (atom_count != total_atoms) {
        fmt::print(stderr, "Warning: expected {} atoms in XYZ frame, read {}\n", total_atoms, atom_count);
    }
}

void System::load_lammpstrj(std::ifstream& file) {
    std::string line;
    total_atoms = 0;

    std::vector<std::string> atom_columns;
    bool got_box = false;
    bool triclinic = false;
    std::array<bool, 3> pbc{{true, true, true}};
    double xlo_bound = 0.0, xhi_bound = 0.0;
    double ylo_bound = 0.0, yhi_bound = 0.0;
    double zlo_bound = 0.0, zhi_bound = 0.0;
    double xy = 0.0, xz = 0.0, yz = 0.0;

    while (std::getline(file, line)) {
        const auto tokens = split_by_space(line);
        if (tokens.empty() || tokens[0] != "ITEM:") continue;
        if (tokens.size() < 2) continue;

        if (tokens[1] == "TIMESTEP") {
            std::getline(file, line);
            continue;
        }
        if (tokens[1] == "NUMBER") {
            if (!std::getline(file, line)) break;
            const auto count_tokens = split_by_space(line);
            if (count_tokens.empty()) break;
            total_atoms = std::stoi(count_tokens[0]);
            atoms.reserve(total_atoms);
            bonds.reserve(total_atoms * 3);
            molecules.reserve(total_atoms / 2);
            continue;
        }
        if (tokens[1] == "BOX") {
            if (tokens.size() < 4 || tokens[2] != "BOUNDS") continue;
            if (tokens[3] == "general") {
                fmt::print(stderr,
                    "Error: general triclinic LAMMPS dump boxes are not supported; "
                    "write a restricted triclinic dump instead\n");
                return;
            }

            triclinic = std::find(tokens.begin(), tokens.end(), "xy") != tokens.end();
            const int boundary_start = 3 + (triclinic ? 3 : 0);
            for (int i = 0; i < 3; ++i) {
                const int pos = boundary_start + i;
                if (pos < static_cast<int>(tokens.size()) && !tokens[pos].empty()) {
                    const char flag = static_cast<char>(std::tolower(static_cast<unsigned char>(tokens[pos][0])));
                    if (flag == 'p') pbc[i] = true;
                    else if (flag == 'f') pbc[i] = false;
                }
            }

            for (int dim = 0; dim < 3; ++dim) {
                if (!std::getline(file, line)) {
                    fmt::print(stderr, "Error: truncated LAMMPS dump BOX BOUNDS\n");
                    return;
                }
                const auto values = split_by_space(line);
                if (values.size() < (triclinic ? 3u : 2u)) {
                    fmt::print(stderr, "Error: invalid LAMMPS dump BOX BOUNDS line\n");
                    return;
                }
                if (dim == 0) {
                    xlo_bound = std::stod(values[0]);
                    xhi_bound = std::stod(values[1]);
                    if (triclinic) xy = std::stod(values[2]);
                }
                else if (dim == 1) {
                    ylo_bound = std::stod(values[0]);
                    yhi_bound = std::stod(values[1]);
                    if (triclinic) xz = std::stod(values[2]);
                }
                else {
                    zlo_bound = std::stod(values[0]);
                    zhi_bound = std::stod(values[1]);
                    if (triclinic) yz = std::stod(values[2]);
                }
            }
            got_box = true;
            continue;
        }
        if (tokens[1] == "ATOMS") {
            atom_columns.assign(tokens.begin() + 2, tokens.end());
            break;
        }
    }

    if (total_atoms <= 0 || atom_columns.empty()) return;
    if (!got_box) {
        fmt::print(stderr, "Error: LAMMPS dump frame has no BOX BOUNDS\n");
        return;
    }

    const double xlo = xlo_bound - std::min({0.0, xy, xz, xy + xz});
    const double xhi = xhi_bound - std::max({0.0, xy, xz, xy + xz});
    const double ylo = ylo_bound - std::min({0.0, yz});
    const double yhi = yhi_bound - std::max({0.0, yz});
    const double zlo = zlo_bound;
    const double zhi = zhi_bound;

    Cell frame_cell;
    frame_cell.set_lammps_restricted(xhi - xlo, yhi - ylo, zhi - zlo, xy, xz, yz, {xlo, ylo, zlo}, pbc);
    cell = frame_cell;
    has_boundaries = cell.valid && cell.any_periodic();
    if (!cell.valid) {
        fmt::print(stderr, "Error: singular LAMMPS dump box\n");
        return;
    }

    auto column_index = [&](const std::string& name) -> int {
        auto found = std::find(atom_columns.begin(), atom_columns.end(), name);
        if (found == atom_columns.end()) return -1;
        return static_cast<int>(found - atom_columns.begin());
    };

    const int id_column = column_index("id");
    const int type_column = column_index("type");
    const int element_column = column_index("element");
    if (type_column < 0 && element_column < 0) {
        fmt::print(stderr, "Error: LAMMPS dump ATOMS has neither type nor element column\n");
        return;
    }

    enum class CoordMode { CARTESIAN, CARTESIAN_UNWRAPPED, SCALED, SCALED_UNWRAPPED };
    CoordMode coord_mode = CoordMode::CARTESIAN;
    int x_column = -1, y_column = -1, z_column = -1;

    const int xsu = column_index("xsu");
    const int ysu = column_index("ysu");
    const int zsu = column_index("zsu");
    const int xs = column_index("xs");
    const int ys = column_index("ys");
    const int zs = column_index("zs");
    const int xu = column_index("xu");
    const int yu = column_index("yu");
    const int zu = column_index("zu");
    const int x = column_index("x");
    const int y = column_index("y");
    const int z = column_index("z");

    if (xsu >= 0 && ysu >= 0 && zsu >= 0) {
        coord_mode = CoordMode::SCALED_UNWRAPPED;
        x_column = xsu;
        y_column = ysu;
        z_column = zsu;
    }
    else if (xs >= 0 && ys >= 0 && zs >= 0) {
        coord_mode = CoordMode::SCALED;
        x_column = xs;
        y_column = ys;
        z_column = zs;
    }
    else if (xu >= 0 && yu >= 0 && zu >= 0) {
        coord_mode = CoordMode::CARTESIAN_UNWRAPPED;
        x_column = xu;
        y_column = yu;
        z_column = zu;
    }
    else if (x >= 0 && y >= 0 && z >= 0) {
        coord_mode = CoordMode::CARTESIAN;
        x_column = x;
        y_column = y;
        z_column = z;
    }
    else {
        fmt::print(stderr, "Error: LAMMPS dump ATOMS has no complete x/y/z column set\n");
        return;
    }

    int required_columns = 0;
    for (int index : {id_column, type_column, element_column, x_column, y_column, z_column}) {
        required_columns = std::max(required_columns, index + 1);
    }

    auto add_type = [&](const std::string& name) -> int {
        auto found = type_stoi.find(name);
        if (found != type_stoi.end()) return found->second;
        const int type_id = System::register_type(name);
        type_itos[type_id] = name;
        type_stoi[name] = type_id;
        total_types = std::max(total_types, type_id);
        return type_id;
    };

    int atoms_count = 0;
    while (atoms_count < total_atoms) {
        if (!std::getline(file, line)) break;
        const auto tokens = split_by_space(line);
        if (tokens.empty()) continue;
        if (static_cast<int>(tokens.size()) < required_columns) {
            fmt::print(stderr, "Error: LAMMPS dump atom line has {} columns, expected at least {}\n",
                tokens.size(), required_columns);
            return;
        }

        const int atom_id = (id_column >= 0) ? std::stoi(tokens[id_column]) : atoms_count + 1;
        int type_id = 0;
        std::string type_name;

        if (type_column >= 0) {
            const std::string raw_type = tokens[type_column];
            if (is_integer_token(raw_type)) {
                type_id = std::stoi(raw_type);
                auto found = type_itos.find(type_id);
                if (found == type_itos.end()) {
                    if (element_column >= 0) {
                        type_name = tokens[element_column];
                        type_id = add_type(type_name);
                    }
                    else {
                        fmt::print(stderr,
                            "Error: LAMMPS dump type {} has no element mapping; provide --types/-t\n", type_id);
                        return;
                    }
                }
                else {
                    type_name = found->second;
                }
            }
            else {
                type_name = raw_type;
                type_id = add_type(type_name);
            }
        }
        else {
            type_name = tokens[element_column];
            type_id = add_type(type_name);
        }

        const double c0 = std::stod(tokens[x_column]);
        const double c1 = std::stod(tokens[y_column]);
        const double c2 = std::stod(tokens[z_column]);
        std::array<double, 3> relative{{0.0, 0.0, 0.0}};

        if (coord_mode == CoordMode::SCALED || coord_mode == CoordMode::SCALED_UNWRAPPED) {
            std::array<double, 3> fractional{{c0, c1, c2}};
            if (coord_mode == CoordMode::SCALED) {
                for (int i = 0; i < 3; ++i) {
                    if (cell.periodic[i]) fractional[i] -= std::floor(fractional[i]);
                }
            }
            relative = cell.cartesian(fractional);
        }
        else {
            relative = {c0 - cell.origin[0], c1 - cell.origin[1], c2 - cell.origin[2]};
            if (coord_mode == CoordMode::CARTESIAN) {
                relative = cell.wrap_relative(relative);
            }
        }

        Atom* atom = new Atom(atom_id, type_id,
            {static_cast<float>(relative[0]), static_cast<float>(relative[1]), static_cast<float>(relative[2])},
            type_name);
        atoms.push_back(atom);
        id_atom_map[atom_id] = atom;
        ++atoms_count;
    }

    if (atoms_count != total_atoms) {
        fmt::print(stderr, "Warning: expected {} atoms in LAMMPS dump frame, read {}\n", total_atoms, atoms_count);
    }
}

/**
 * @brief Set atom types and initialize type-related mappings and valences.
 * @param type_names Vector of atom type names.
 */
void Molecule::update_topology() {
    types_nums.clear();
    std::string mol_formula;
    for (auto& mol_atom : mol_atoms) {
        types_nums[mol_atom->type_name]++;
    }
    for (auto& pair : types_nums) {
        mol_formula += fmt::format("{}{}", pair.first, pair.second);
    }
    formula = rename_formula(mol_formula, ELEMENT_DISPLAY_ORDER);

    // >200 atoms: formula-level identity, no SMILES (kerogen/polymer regime).
    if (mol_atoms.size() > reax_smiles::kMaxSmilesAtoms) {
        smiles.clear();
        structure_id = reax_smiles::formula_id(formula);
        hash = reax_string_hash(formula);
        for (auto& atom : mol_atoms) atom->belong_molecule = this;
        return;
    }

    std::vector<Atom*> atoms(mol_atoms.begin(), mol_atoms.end());
    std::sort(atoms.begin(), atoms.end(),
        [](const Atom* a, const Atom* b) { return a->id < b->id; });

    std::map<const Atom*, int> atom_index;
    for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
        atom_index[atoms[i]] = i;
    }

    std::map<const Atom*, int> degree;
    for (const Bond* bond : mol_bonds) {
        degree[bond->atom_i]++;
        degree[bond->atom_j]++;
    }

    std::vector<reax_smiles::Atom> smiles_atoms;
    smiles_atoms.reserve(atoms.size());
    for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
        const std::string& element = atoms[i]->type_name;
        int valence = degree[atoms[i]];
        // Formal charges are only used to keep standard-valence SMILES
        // parseable; they are not a charge-equilibration model. ReaxFF
        // ammonium (N valence 4) and oxonium (O valence 3) are common.
        int charge = 0;
        if (element == "N" && valence > 3) charge = valence - 3;
        else if (element == "O" && valence > 2) charge = valence - 2;
        smiles_atoms.push_back({i, element, charge});
    }

    std::vector<reax_smiles::Bond> smiles_bonds;
    smiles_bonds.reserve(mol_bonds.size());
    for (const Bond* bond : mol_bonds) {
        auto i = atom_index.find(bond->atom_i);
        auto j = atom_index.find(bond->atom_j);
        if (i == atom_index.end() || j == atom_index.end()) continue;
        smiles_bonds.push_back({i->second, j->second});
    }
    std::sort(smiles_bonds.begin(), smiles_bonds.end(), [](const auto& a, const auto& b) {
        if (a.i != b.i) return a.i < b.i;
        return a.j < b.j;
    });

    reax_smiles::Result result = reax_smiles::canonical_smiles(smiles_atoms, smiles_bonds);
    if (result.fallback || result.smiles.empty()) {
        smiles.clear();
        structure_id = reax_smiles::formula_id(formula);
    } else {
        smiles = result.smiles;
        structure_id = reax_smiles::smiles_id(smiles);
    }

    hash = reax_string_hash(formula);  // legacy formula hash, kept for now
    for (auto& atom : mol_atoms) {
        atom->belong_molecule = this;
    }
}

void System::reset_type_registry() {
    global_type_stoi.clear();
    global_type_itos.clear();
    int type_id = 0;
    for (const auto& name : INPUT_ELEMENT_TYPES) {
        if (name.empty() || global_type_stoi.find(name) != global_type_stoi.end()) continue;
        ++type_id;
        global_type_stoi[name] = type_id;
        global_type_itos[type_id] = name;
    }
    global_type_registry_initialized = true;
}

int System::register_type(const std::string& name) {
    if (!global_type_registry_initialized) reset_type_registry();
    auto found = global_type_stoi.find(name);
    if (found != global_type_stoi.end()) return found->second;

    int next_id = 1;
    if (!global_type_itos.empty()) next_id = global_type_itos.rbegin()->first + 1;
    global_type_stoi[name] = next_id;
    global_type_itos[next_id] = name;
    return next_id;
}

void System::set_types() {
    if (!global_type_registry_initialized) reset_type_registry();

    type_stoi = global_type_stoi;
    type_itos = global_type_itos;
    total_types = type_itos.empty() ? 0 : type_itos.rbegin()->first;

    for (const auto& pair : type_itos) {
        const std::string& name = pair.second;
        // Initialize ELEMENT_MAX_VALENCIES for each type.
        if (ELEMENT_MAX_VALENCIES.find(name) == ELEMENT_MAX_VALENCIES.end()) {
            ELEMENT_MAX_VALENCIES[name] = 4;  // Default max valence
        }
    }
}

/**
 * @brief Print system summary for the current frame, including atom, bond, molecule, and ring counts.
 */
void System::finish() {
    std::string header = fmt::format("{:<8} {:<8} {:<8} {:<8}", "Frame", "Atoms", "Bonds", "Mols");
    std::string formatting = "{:<8} {:<8} {:<8} {:<8} {}";

    std::string ring_str = "";

    if (!FLAG_NO_RINGS) {
        for (int i = MIN_RING_SIZE; i <= MAX_RING_SIZE; i++) {
            header += fmt::format(" R{:<5}", i);

            if (ring_counts.find(i) != ring_counts.end()) {
                ring_str += fmt::format(" {:<6}", ring_counts[i]);
            }
            else {
                ring_str += fmt::format(" {:<6}", 0);
            }
        }
    }

    bool should_print_frame = frame_id <= 10 || (frame_id <= 100 && frame_id % 10 == 0) || frame_id % 100 == 0;

    if (frame_id <= 1) {
        fmt::print(header);
        fmt::print("\n");
    }

    if (should_print_frame) {
        fmt::print(formatting, frame_id, atoms.size(), bonds.size(), molecules.size(), ring_str);
        fmt::print("\n");
    }

    static int has_warned_pbc = 5;
    if (!this->has_boundaries && has_warned_pbc > 0) {
        fmt::print("**** WARNING **** : No PBC info, cost much more time and results may be inaccurate!\n");
        has_warned_pbc--;
    }
    else if (!this->has_boundaries && has_warned_pbc == 0) {
        fmt::print("**** WARNING **** : No more PBC warnings, you'd better know what you are doing!\n");
        has_warned_pbc--;
    }

    static int has_warned_atom_lost = 5;
    if (prev_sys != nullptr) {
        if (prev_sys->atoms.size() != atoms.size() && has_warned_atom_lost > 0) {
            fmt::print("**** WARNING **** : Detected ATOM LOST, reactions may be inaccurate!\n");
            has_warned_atom_lost--;
        }
        else if (prev_sys->atoms.size() != atoms.size() && has_warned_atom_lost == 0) {
            fmt::print("**** WARNING **** : No more ATOM LOST warnings, you'd better know what you are doing!\n");
            has_warned_atom_lost--;
        }
    }
}

/**
 * @brief Naive neighbor search for systems without periodic boundaries.
 *        Populates each atom's neighbor list based on a cutoff radius.
 */
void System::search_neigh_naive() {
    float radius_sq = neigh_radius * neigh_radius;
    float dist_sq;
    // Naive search for test

    // Define a function pointer type for distance calculation function
    for (auto& curr_atom : atoms) {
        for (auto& other_atom : atoms) {
            if (curr_atom == other_atom) continue;

            dist_sq = cell.valid ? cell.distance_sq(curr_atom->coord, other_atom->coord)
                                 : distance_sq(curr_atom->coord, other_atom->coord);
            if (dist_sq <= radius_sq) {
                curr_atom->neighs.insert(other_atom);
            }
            if (curr_atom->neighs.size() >= MAX_NEIGH) {
                continue;
            }
        }
    }
}

/**
 * @brief Use cell list algorithm to search for neighbors within a cutoff radius.
 */
void System::search_neigh_cell_list() {
    Cell_list cell_list(atoms, neigh_radius, cell, MAX_NEIGH);
    for (auto& atom : atoms) {
        cell_list.search_neighbors(atom);
    }
}

/**
 * @brief Selects the appropriate neighbor search method based on boundary conditions.
 */
void System::search_neigh() {
    if (cell.valid && cell.periodic[0] && cell.periodic[1] && cell.periodic[2]) {
        search_neigh_cell_list();
    }
    else {
        search_neigh_naive();
    }
}

/**
 * @brief Build bonds between atoms based on distance and type-dependent radii.
 *        Ensures valence limits and updates bond/atom connectivity.
 */
void System::build_bonds_by_radius() {
    // prepare types and radius.

    // compute.
    float bond_r;
    float bond_sq;
    float dist_sq;
    float relative_sq;

    Atom* tmp_neigh;
    std::vector<std::pair<Atom*, float>> candidate_id_relative_sq;
    candidate_id_relative_sq.reserve(MAX_NEIGH);

    std::pair<int, int> type_idx;

    for (auto& atom : atoms) {
        // Skip X atoms and atoms that reached max valence
        if (type_itos[atom->type_id] == "X" || atom->bonded_atoms.size() >= atom->max_valence) continue;

        for (auto& neigh : atom->neighs) {
            type_idx = { atom->type_id, neigh->type_id };
            bond_sq = bond_radius_sq[type_idx];

            if (bond_sq <= 0) continue;

            if (cell.valid) {
                dist_sq = cell.distance_sq(atom->coord, neigh->coord);
            }
            else {
                dist_sq = distance_sq(atom->coord, neigh->coord);
            }

            if (dist_sq > bond_sq) continue;

            relative_sq = dist_sq / bond_sq;

            if (relative_sq < 1.0f) {
                candidate_id_relative_sq.emplace_back(std::pair(neigh, relative_sq));
            }
        }

        std::sort(candidate_id_relative_sq.begin(), candidate_id_relative_sq.end(),
            [](const auto& a, const auto& b) { return a.second < b.second; });

        for (size_t tmp_id = 0; tmp_id < candidate_id_relative_sq.size(); tmp_id++) {
            if (atom->bonded_atoms.size() >= atom->max_valence) break;

            tmp_neigh = candidate_id_relative_sq[tmp_id].first;

            if (atom->bonded_atoms.find(tmp_neigh) != atom->bonded_atoms.end()) continue;

            if (tmp_neigh->bonded_atoms.size() >= tmp_neigh->max_valence) continue;

            Bond* bond = new Bond(atom, tmp_neigh);
            bond_type_counts[std::pair(atom->type_id, tmp_neigh->type_id)]++;

            bonds.emplace_back(bond);
            atom->bonds.insert(bond);
            tmp_neigh->bonds.insert(bond);
            atom->bonded_atoms.insert(tmp_neigh);
            tmp_neigh->bonded_atoms.insert(atom);
        }

        candidate_id_relative_sq.clear();
    }
}

/**
 * @brief Build molecules by grouping bonded atoms and maximizing bond orders within valence constraints.
 */
void System::build_molecules() {
    std::set<Atom*> visited;
    int new_mol_id = 0;

    for (auto& atom : atoms) {
        if (visited.find(atom) == visited.end())  // If not visited
        {
            new_mol_id++;
            Molecule* new_mol = new Molecule(new_mol_id);
            dfs(atom, visited, new_mol);
            molecules.push_back(new_mol);
        }
    }

    // Weighing bonds
    // Maximize total weight of edges (bond orders)
    // Constraint: each atom's valence can not exceeded.
    std::map<Atom*, int> atom_valence_residual;
    for (auto& atom : atoms) {
        atom_valence_residual[atom] = atom->max_valence - atom->bonded_atoms.size();
    }

    int current_weight_gain = 0;

    for (auto& molecule : molecules) {
        // Skip single atom molecules and molecules containing X atoms
        if (molecule->mol_atoms.size() == 1) continue;

        // Greedy algorithm.
        while (true) {
            current_weight_gain = 0;

            for (auto& bond : molecule->mol_bonds) {
                if (atom_valence_residual[bond->atom_i] > 0 && atom_valence_residual[bond->atom_j] > 0) {
                    current_weight_gain += 2;
                    bond->order++;
                    atom_valence_residual[bond->atom_i]--;
                    atom_valence_residual[bond->atom_j]--;

                    // We do not consider 4+ bonds, they merely exist.
                    if (bond->order >= 3) continue;
                }
            }

            if (current_weight_gain == 0) {
                break;
            }
        }
    }

    atom_valence_residual.clear();

    for (auto& molecule : molecules) {
        molecule->update_topology();
        id_mol_map[molecule->id] = molecule;
    }

    visited.clear();
}

/**
 * @brief Depth-first search to collect all atoms and bonds in a molecule starting from a given atom.
 * @param atom Starting atom for DFS.
 * @param visited Set of already visited atoms.
 * @param curr_mol Molecule object being constructed.
 */
void System::dfs(Atom* atom, std::set<Atom*>& visited, Molecule* curr_mol) {
    visited.insert(atom);
    curr_mol->mol_atoms.insert(atom);
    curr_mol->atom_ids.insert(atom->id);

    for (auto& bond : atom->bonds) {
        curr_mol->mol_bonds.insert(bond);
        // One of atom_i and atom_j is the atom itself.
        if (visited.find(bond->atom_i) == visited.end()) {
            dfs(bond->atom_i, visited, curr_mol);
            continue;
        }
        if (visited.find(bond->atom_j) == visited.end()) {
            dfs(bond->atom_j, visited, curr_mol);
        }
    }
}

/**
 * @brief Dump the current system as a LAMMPS data file for testing.
 *
 * Molecule-ID column: by default each independent molecule gets an id sorted by
 * atom count (largest first).  With --mark-ring-atoms the column is instead the
 * atom's smallest ring size, matching the legacy "use molid for ring marking"
 * behaviour.  Angles and dihedrals are written as topology-only sections.
 */
void System::dump_lammps_data() {
    if ((frame_id - 1) % DUMP_STEPS != 0) return;

    const std::string file_basename = fmt::format("frame_{}.data", frame_id);
    FILE* file = create_file(file_basename);
    if (file == nullptr) return;

    std::set<std::array<int, 3>> angles;
    std::set<std::array<int, 4>> dihedrals;

    auto neighboring_ids = [](Atom* atom, int exclude_id) {
        std::vector<int> ids;
        ids.reserve(atom->bonds.size());
        for (Bond* bond : atom->bonds) {
            Atom* other = (bond->atom_i == atom) ? bond->atom_j : bond->atom_i;
            if (other->id != exclude_id) ids.push_back(other->id);
        }
        std::sort(ids.begin(), ids.end());
        return ids;
    };

    for (Atom* atom : atoms) {
        const auto ids = neighboring_ids(atom, -1);
        for (std::size_t i = 0; i < ids.size(); ++i) {
            for (std::size_t j = i + 1; j < ids.size(); ++j) {
                angles.insert({ids[i], atom->id, ids[j]});
            }
        }
    }
    for (Bond* bond : bonds) {
        const auto left = neighboring_ids(bond->atom_i, bond->atom_j->id);
        const auto right = neighboring_ids(bond->atom_j, bond->atom_i->id);
        for (int a : left) {
            for (int b : right) {
                std::array<int, 4> forward{{a, bond->atom_i->id, bond->atom_j->id, b}};
                std::array<int, 4> reverse{{b, bond->atom_j->id, bond->atom_i->id, a}};
                dihedrals.insert(std::min(forward, reverse));
            }
        }
    }

    const auto dump_ids = molecule_dump_ids();
    Cell data_cell;
    std::array<std::array<double, 3>, 3> rotation{{{{1.0, 0.0, 0.0}}, {{0.0, 1.0, 0.0}}, {{0.0, 0.0, 1.0}}}};
    bool have_cell = false;
    std::array<double, 3> shift{{0.0, 0.0, 0.0}};

    if (cell.valid) {
        have_cell = cell.to_restricted(data_cell, rotation);
        if (!have_cell) {
            fmt::print(stderr, "Warning: could not convert cell to LAMMPS restricted form; using bounding box\n");
        }
    }

    if (!have_cell) {
        std::array<double, 3> lower{{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                                     std::numeric_limits<double>::max()}};
        std::array<double, 3> upper{{-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
                                     -std::numeric_limits<double>::max()}};
        for (Atom* atom : atoms) {
            for (int i = 0; i < 3; ++i) {
                lower[i] = std::min(lower[i], static_cast<double>(atom->coord[i]));
                upper[i] = std::max(upper[i], static_cast<double>(atom->coord[i]));
            }
        }
        if (atoms.empty()) {
            lower = {0.0, 0.0, 0.0};
            upper = {1.0, 1.0, 1.0};
        }
        const double margin = 0.5;
        std::array<double, 3> lengths{{0.0, 0.0, 0.0}};
        for (int i = 0; i < 3; ++i) {
            lengths[i] = std::max(upper[i] - lower[i] + 2.0 * margin, 2.0 * margin);
            shift[i] = -lower[i] + margin;
        }
        data_cell.set_lattice_diagonal(lengths, {0.0, 0.0, 0.0}, {false, false, false});
    }

    auto data_coord = [&](Atom* atom) {
        std::array<double, 3> coord{{static_cast<double>(atom->coord[0]), static_cast<double>(atom->coord[1]),
                                     static_cast<double>(atom->coord[2])}};
        if (have_cell) {
            const auto H = rotation;
            return std::array<double, 3>{{H[0][0] * coord[0] + H[0][1] * coord[1] + H[0][2] * coord[2],
                                          H[1][0] * coord[0] + H[1][1] * coord[1] + H[1][2] * coord[2],
                                          H[2][0] * coord[0] + H[2][1] * coord[1] + H[2][2] * coord[2]}};
        }
        return std::array<double, 3>{{coord[0] + shift[0], coord[1] + shift[1], coord[2] + shift[2]}};
    };

    int ntypes = total_types;
    for (const auto& pair : type_itos) ntypes = std::max(ntypes, pair.first);
    ntypes = std::max(ntypes, 1);

    fmt::print(file, "# LAMMPS data file written by reax_tools\n\n");
    fmt::print(file, "{} atoms\n", atoms.size());
    fmt::print(file, "{} bonds\n", bonds.size());
    fmt::print(file, "{} angles\n", angles.size());
    fmt::print(file, "{} dihedrals\n", dihedrals.size());
    fmt::print(file, "{} atom types\n", ntypes);
    fmt::print(file, "1 bond types\n");
    fmt::print(file, "1 angle types\n");
    fmt::print(file, "1 dihedral types\n\n");

    const double lx = data_cell.a[0];
    const double ly = data_cell.b[1];
    const double lz = data_cell.c[2];
    const double xy = data_cell.b[0];
    const double xz = data_cell.c[0];
    const double yz = data_cell.c[1];
    if (data_cell.is_restricted_orthogonal()) {
        fmt::print(file, "{} {} xlo xhi\n", 0.0, lx);
        fmt::print(file, "{} {} ylo yhi\n", 0.0, ly);
        fmt::print(file, "{} {} zlo zhi\n", 0.0, lz);
    }
    else {
        fmt::print(file, "{} {} {} xlo xhi xy\n", 0.0, lx, xy);
        fmt::print(file, "{} {} {} ylo yhi xz\n", 0.0, ly, xz);
        fmt::print(file, "{} {} {} zlo zhi yz\n", 0.0, lz, yz);
    }

    fmt::print(file, "\nMasses\n\n");
    for (int type_id = 1; type_id <= ntypes; ++type_id) {
        std::string symbol;
        auto found = type_itos.find(type_id);
        if (found != type_itos.end()) symbol = found->second;
        auto mass_it = ELEMENT_MASS.find(symbol);
        const double mass = (mass_it == ELEMENT_MASS.end()) ? 1.0 : mass_it->second;
        fmt::print(file, "{} {:.6f} # {}\n", type_id, mass, symbol.empty() ? "X" : symbol);
    }

    fmt::print(file, "\nAtoms # full\n\n");
    if (FLAG_MARK_RING_ATOMS) {
        std::vector<Atom*> sorted_atoms = atoms;
        std::sort(sorted_atoms.begin(), sorted_atoms.end(),
            [](const Atom* a, const Atom* b) { return a->id < b->id; });
        for (Atom* atom : sorted_atoms) {
            const auto coord = data_coord(atom);
            fmt::print(file, "{} {} {} 0.0 {:>.6f} {:>.6f} {:>.6f} 0 0 0\n", atom->id, atom->ringsize,
                atom->type_id, coord[0], coord[1], coord[2]);
        }
    }
    else {
        for (Molecule* molecule : molecules_by_size()) {
            const int molid = dump_ids.at(molecule);
            std::vector<Atom*> sorted_atoms(molecule->mol_atoms.begin(), molecule->mol_atoms.end());
            std::sort(sorted_atoms.begin(), sorted_atoms.end(),
                [](const Atom* a, const Atom* b) { return a->id < b->id; });
            for (Atom* atom : sorted_atoms) {
                const auto coord = data_coord(atom);
                fmt::print(file, "{} {} {} 0.0 {:>.6f} {:>.6f} {:>.6f} 0 0 0\n", atom->id, molid,
                    atom->type_id, coord[0], coord[1], coord[2]);
            }
        }
    }

    fmt::print(file, "\nBonds\n\n");
    for (std::size_t i = 0; i < bonds.size(); ++i) {
        fmt::print(file, "{} 1 {} {}\n", i + 1, bonds[i]->atom_i->id, bonds[i]->atom_j->id);
    }

    if (!angles.empty()) {
        fmt::print(file, "\nAngles\n\n");
        int angle_id = 0;
        for (const auto& angle : angles) {
            fmt::print(file, "{} 1 {} {} {}\n", ++angle_id, angle[0], angle[1], angle[2]);
        }
    }

    if (!dihedrals.empty()) {
        fmt::print(file, "\nDihedrals\n\n");
        int dihedral_id = 0;
        for (const auto& dihedral : dihedrals) {
            fmt::print(file, "{} 1 {} {} {} {}\n", ++dihedral_id, dihedral[0], dihedral[1], dihedral[2],
                dihedral[3]);
        }
    }

    fmt::print(file, "\n");
    fclose(file);
}

std::vector<Molecule*> System::molecules_by_size() const {
    std::vector<Molecule*> sorted = molecules;
    auto min_atom_id = [](const Molecule* molecule) {
        int result = std::numeric_limits<int>::max();
        for (int id : molecule->atom_ids) result = std::min(result, id);
        return result;
    };
    std::sort(sorted.begin(), sorted.end(), [&](const Molecule* a, const Molecule* b) {
        if (a->mol_atoms.size() != b->mol_atoms.size()) return a->mol_atoms.size() > b->mol_atoms.size();
        if (a->formula != b->formula) return a->formula < b->formula;
        return min_atom_id(a) < min_atom_id(b);
    });
    return sorted;
}

std::map<const Molecule*, int> System::molecule_dump_ids() const {
    std::map<const Molecule*, int> ids;
    int next_id = 1;
    for (const Molecule* molecule : molecules_by_size()) {
        ids[molecule] = next_id++;
    }
    return ids;
}

void System::dump_xyz_marked(FILE* file, const std::vector<std::string>& element_columns) {
    if (file == nullptr) return;
    if ((frame_id - 1) % DUMP_XYZ_STEPS != 0) return;

    const auto dump_ids = molecule_dump_ids();
    std::map<const Molecule*, double> molecule_mass;
    for (const Molecule* molecule : molecules) {
        double mass = 0.0;
        for (const Atom* atom : molecule->mol_atoms) {
            auto found = ELEMENT_MASS.find(atom->type_name);
            mass += (found == ELEMENT_MASS.end()) ? 0.0 : found->second;
        }
        molecule_mass[molecule] = mass;
    }

    Cell output_cell;
    std::array<double, 3> shift{{0.0, 0.0, 0.0}};
    bool preserve_origin = false;

    if (cell.valid) {
        output_cell = cell;
        preserve_origin = std::fabs(cell.origin[0]) > 1.0e-12 || std::fabs(cell.origin[1]) > 1.0e-12 ||
                          std::fabs(cell.origin[2]) > 1.0e-12;
    }
    else {
        std::array<double, 3> lower{{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                                     std::numeric_limits<double>::max()}};
        std::array<double, 3> upper{{-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
                                     -std::numeric_limits<double>::max()}};
        for (const Atom* atom : atoms) {
            for (int i = 0; i < 3; ++i) {
                const double value = static_cast<double>(atom->coord[i]) + cell.origin[i];
                lower[i] = std::min(lower[i], value);
                upper[i] = std::max(upper[i], value);
            }
        }
        if (atoms.empty()) {
            lower = {0.0, 0.0, 0.0};
            upper = {1.0, 1.0, 1.0};
        }
        const double margin = 0.5;
        std::array<double, 3> lengths{{0.0, 0.0, 0.0}};
        for (int i = 0; i < 3; ++i) {
            lengths[i] = std::max(upper[i] - lower[i] + 2.0 * margin, 2.0 * margin);
            shift[i] = -lower[i] + margin;
        }
        output_cell.set_lattice_diagonal(lengths, {0.0, 0.0, 0.0}, {false, false, false});
    }

    std::string properties = "species:S:1:pos:R:3:id:I:1:molecule:I:1:molsize:I:1:ringsize:I:1:formula:S:1:molmass:R:1";
    for (const auto& element : element_columns) {
        properties += fmt::format(":mol{}:I:1", element);
    }

    fmt::print(file, "{}\n", atoms.size());
    fmt::print(file,
        "Lattice=\"{:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f}\" "
        "Properties={} pbc=\"{} {} {}\"",
        output_cell.a[0], output_cell.a[1], output_cell.a[2], output_cell.b[0], output_cell.b[1],
        output_cell.b[2], output_cell.c[0], output_cell.c[1], output_cell.c[2], properties,
        output_cell.periodic[0] ? "T" : "F", output_cell.periodic[1] ? "T" : "F",
        output_cell.periodic[2] ? "T" : "F");
    if (preserve_origin) {
        fmt::print(file, " Origin=\"{:.8f} {:.8f} {:.8f}\"", cell.origin[0], cell.origin[1], cell.origin[2]);
    }
    fmt::print(file, " frame={}\n", frame_id);

    for (const Atom* atom : atoms) {
        const Molecule* molecule = atom->belong_molecule;
        const int molid = molecule ? dump_ids.at(molecule) : 0;
        const int molsize = molecule ? static_cast<int>(molecule->mol_atoms.size()) : 1;
        const double molmass = molecule ? molecule_mass.at(molecule) : 0.0;
        std::string formula = molecule ? molecule->formula : "unknown";
        if (formula.empty()) formula = "unknown";

        std::array<double, 3> coord{{static_cast<double>(atom->coord[0]), static_cast<double>(atom->coord[1]),
                                     static_cast<double>(atom->coord[2])}};
        if (cell.valid) {
            coord[0] += cell.origin[0];
            coord[1] += cell.origin[1];
            coord[2] += cell.origin[2];
        }
        else {
            coord[0] += shift[0];
            coord[1] += shift[1];
            coord[2] += shift[2];
        }

        fmt::print(file, "{} {:.8f} {:.8f} {:.8f} {} {} {} {} {} {:.6f}", atom->type_name, coord[0], coord[1],
            coord[2], atom->id, molid, molsize, atom->ringsize, formula, molmass);
        for (const auto& element : element_columns) {
            int count = 0;
            if (molecule) {
                auto found = molecule->types_nums.find(element);
                if (found != molecule->types_nums.end()) count = found->second;
            }
            fmt::print(file, " {}", count);
        }
        fmt::print(file, "\n");
    }
}

/**
 * @brief Process the current system: neighbor search, bond building, molecule construction, and ring counting.
 */
void System::process_this() {
    // define bond radius for each type.

    std::map<int, float> atomic_radius;
    float bond_radius = 0.0f;

    for (auto& pair : type_stoi) {
        std::string typ_s = pair.first;
        int typ_i = pair.second;

        // pair.first: string type, pair.second: int type
        if (ELEMENT_ATOMIC_RADII.find(typ_s) != ELEMENT_ATOMIC_RADII.end())
            atomic_radius[typ_i] = ELEMENT_ATOMIC_RADII.at(typ_s);
        else
            atomic_radius[typ_i] = ELEMENT_ATOMIC_RADII.at("X");
    }

    for (int type_i = 1; type_i <= total_types; type_i++) {
        for (int type_j = 1; type_j <= total_types; type_j++) {
            std::pair<int, int> pair_ij = { type_i, type_j };
            if (atomic_radius[type_i] == 0.0f || atomic_radius[type_j] == 0.0f) {
                bond_radius = 0.0f;
            }
            else {
                bond_radius = 0.5f * (atomic_radius[type_i] + atomic_radius[type_j]) * RVDW_FACTOR;
            }
            if (bond_radius > neigh_radius) {
                neigh_radius = bond_radius;
            }
            bond_radius_sq[pair_ij] = bond_radius * bond_radius;
            bond_type_counts[pair_ij] = 0;
        }
    }

    search_neigh();
    build_bonds_by_radius();
    build_molecules();

    if (!FLAG_NO_RINGS) {
        compute_ring_counts();
    }
}

/**
 * @brief Process the current system frame for ReaxFF analysis, including formula extraction and reaction tracking.
 *        Thread-safe import of frame formulas and molecule comparison with previous frame.
 */
void System::process_counters() {
    std::vector<SpeciesRecord> frame_species(molecules.size());
    std::map<std::string, int> bond_count_map;
    std::map<std::string, int> ring_count_map;
    std::map<std::string, int> atom_bonded_num_count_map;

    // Structural species records: id + formula + canonical SMILES.
    for (size_t i = 0; i < molecules.size(); i++) {
        frame_species[i] = {
            molecules[i]->structure_id,
            molecules[i]->formula,
            molecules[i]->smiles,
        };
    }

    // Bonds
    for (const auto& [ij, count] : bond_type_counts) {
        std::string type_a = std::max(type_itos[ij.first], type_itos[ij.second]);
        std::string type_b = std::min(type_itos[ij.first], type_itos[ij.second]);

        std::string key = fmt::format("{}-{}", type_a, type_b);

        if (bond_count_map.find(key) == bond_count_map.end()) {
            bond_count_map[key] = count;
        }
        else {
            bond_count_map[key] += count;
        }
    }

    // Atom bonded nums
    std::map<std::pair<int, int>, int> tmp_atom_bonded_num_count;
    for (const auto& [type_id, type_string] : type_itos) {
        int max_valence = ELEMENT_MAX_VALENCIES[type_string];
        for (size_t valence = 0; valence <= max_valence; valence++) {
            tmp_atom_bonded_num_count[{type_id, valence}] = 0;
        }
    }
    for (const auto& atom : atoms) {
        tmp_atom_bonded_num_count[{atom->type_id, atom->bonded_atoms.size()}]++;
    }

    for (const auto& [id_num, count] : tmp_atom_bonded_num_count) {
        std::string key = fmt::format("{}-{}", type_itos[id_num.first], id_num.second);
        atom_bonded_num_count_map[key] = count;
    }

    // Rings
    for (const auto& [ring_size, count] : ring_counts) {
        std::string key = fmt::format("ring size {}", ring_size);
        ring_count_map[key] = count;
    }

    // Lock the shared species_counter import operations.
    {
        std::lock_guard<std::mutex> lock(species_counter_mutex);
        species_counter->import_frame_species(frame_id, frame_species);
        bond_counter->import_frame_values(frame_id, bond_count_map);
        atom_bonded_num_counter->import_frame_values(frame_id, atom_bonded_num_count_map);
        ring_counter->import_frame_values(frame_id, ring_count_map);
    }
}

inline Atom* System::get_atom_by_id(const int& atom_id) {
    // Use id_atom_map to find the atom by its id. Return nullptr if not found.
    auto it = id_atom_map.find(atom_id);
    if (it != id_atom_map.end()) {
        return it->second;
    }
    else {
        return nullptr;
    }
}

inline Molecule* System::get_molecule_by_id(const int& mol_id) {
    // Use id_atom_map to find the atom by its id. Return nullptr if not found.
    auto it = id_mol_map.find(mol_id);
    if (it != id_mol_map.end()) {
        return it->second;
    }
    else {
        return nullptr;
    }
}

void System::compute_ring_counts() {
    // Initialize ring counts for sizes 3 to MAX_RING_SIZE
    // MAX_RING_SIZE in defines.h / defines.cpp
    for (int i = MIN_RING_SIZE; i <= MAX_RING_SIZE; i++) {
        ring_counts[i] = 0;
    }

    std::unordered_set<Atom*> visited;
    std::unordered_set<std::unordered_set<Atom*>*> current_rings;
    std::vector<Atom*> current_path;

    // Simple ring detection algorithm.
    for (auto& molecule : molecules) {
        visited.clear();
        current_rings.clear();
        current_path.clear();

        // At least 3 atoms to form rings.
        if (molecule->mol_atoms.size() < 3) continue;

        // Search from every atom is mandatory for polycyclics like Naphthalene.
        // Actually this is not very time-consuming.
        for (auto& start_atom : molecule->mol_atoms) {
            find_rings_from_atom(start_atom, start_atom, 0, visited, current_rings, current_path);
        }

        // Set ring-size property for all atoms in rings.  An atom in fused
        // rings keeps the smallest ring size containing it.
        for (auto& ring : current_rings) {
            const int ring_size = static_cast<int>(ring->size());
            ring_counts[ring_size]++;
            for (auto& atom : *ring) {
                atom->on_ring = true;
                if (atom->ringsize == 0 || ring_size < atom->ringsize) atom->ringsize = ring_size;
            }
            delete ring;  // Only count matters.
        }
    }

    visited.clear();
    current_rings.clear();
    current_path.clear();
}

/**
 * @brief Computes the number of rings of different sizes in the system using DFS-based ring detection.
 *        Updates the ring_counts member variable.
 */
void System::find_rings_from_atom(Atom* current, Atom* start, int depth, std::unordered_set<Atom*>& visited,
    std::unordered_set<std::unordered_set<Atom*>*>& current_rings,
    std::vector<Atom*>& current_path) {
    // Move forward
    // Add current atom to path
    // depth start from 0 and RING_SIZE start from 1.
    current_path.push_back(current);
    visited.insert(current);
    std::unordered_set<std::unordered_set<Atom*>*> to_erase;
    std::unordered_set<Atom*>* this_ring;

    for (Atom* bonded_atom : current->bonded_atoms) {
        // If ring closed and size = MIN - MAX
        if (bonded_atom == start && depth >= MIN_RING_SIZE - 1 && depth <= MAX_RING_SIZE - 1) {
            this_ring = new std::unordered_set<Atom*>();  // Initialize the pointer
            for (Atom* atom : current_path) {
                this_ring->insert(atom);
            }

            // Determine how to handle this new ring.
            bool to_insert = true;
            if (current_rings.size() > 0) {
                for (auto& other_ring : current_rings) {
                    bool is_subset_of_other = true;
                    bool is_subset_of_this = true;

                    // When this ring has an unique atom, this ring is not a subset of the other one.
                    for (auto& this_atom : *this_ring) {
                        if (other_ring->find(this_atom) == other_ring->end()) {
                            is_subset_of_other = false;
                            break;
                        }
                    }

                    // Vice versa
                    for (auto& other_atom : *other_ring) {
                        if (this_ring->find(other_atom) == this_ring->end()) {
                            is_subset_of_this = false;
                            break;
                        }
                    }

                    // If ring is a proper subset of another, we drop the bigger one.
                    if (is_subset_of_other && !is_subset_of_this) {  // The other contains this one
                        to_erase.insert(other_ring);
                        // Don't change to_insert here
                    }
                    else if (is_subset_of_this && !is_subset_of_other) {  // This one contains the other
                        to_insert = false;
                        break;                                             // No need to check other rings
                    }
                    else if (is_subset_of_this && is_subset_of_other) {  // They are the same
                        to_insert = false;
                        break;  // No need to check other rings
                    }
                    // else: None of them is subset, may have intersection, but that is what we want.
                    // to_insert stays true
                }
            }

            for (auto& erase_ring : to_erase) {
                current_rings.erase(erase_ring);
                delete erase_ring;
            }
            to_erase.clear();

            if (to_insert) {
                current_rings.insert(this_ring);
            }
            else {
                // Free memory if we're not inserting this ring
                delete this_ring;
            }
        }
        else if (visited.find(bonded_atom) == visited.end() && depth < MAX_RING_SIZE - 1) {
            // recursive, if not visited
            find_rings_from_atom(bonded_atom, start, depth + 1, visited, current_rings, current_path);
        }
    }

    // Backtrack: remove current atom from path and visited set
    current_path.pop_back();
    visited.erase(current);
}
