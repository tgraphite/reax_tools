#include "fmt/format.h"
#include "string_tools.h"
#include "universe.h"

void System::load_xyz(std::ifstream& file) {
    std::string line;
    std::vector<std::string> tokens;
    int atom_id = 0;
    int atom_count = 0;

    float xlo, xhi, ylo, yhi, zlo, zhi, lx, ly, lz;
    float x, y, z;

    std::string given_type;
    int type_int;
    std::string type_str;
    std::string tmp_str;

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file");
    }

    while (getline(file, line)) {
        tokens = split_by_space(line);

        if ((tokens.size() == 0) or (tokens[0] == "#")) {
            continue;
        }
        // Skip time or energy information lines
        else if (std::find(tokens.begin(), tokens.end(), "time") != tokens.end() &&
                 std::find(tokens.begin(), tokens.end(), "=") != tokens.end() &&
                 std::find(tokens.begin(), tokens.end(), "E") != tokens.end()) {
            this->has_boundaries = false;
            continue;
        }
        // Atom numbers line
        else if (tokens.size() == 1 && can_convert_to_int(tokens[0])) {
            iatoms = std::stoi(tokens[0]);
            atoms.reserve(iatoms);
            bonds.reserve(iatoms * 3);
            molecules.reserve(iatoms / 2);
            continue;
        } else if (starts_with(tokens[0], "Lattice")) {
            // Lattice="25.774014 0.0 0.0 0.0 65.670043 0.0 0.0 0.0 85.0" Origin="-0.058304 0.0 30.0"...

            // Lattice="10.000
            tmp_str = tokens[0];
            lx = std::stof(split(tmp_str, "\"")[1]);

            // 10.000
            ly = std::stof(tokens[4]);

            // 10.000"
            tmp_str = tokens[8];
            lz = std::stof(split(tmp_str, "\"")[0]);

            // Origin="0.000
            tmp_str = tokens[9];
            xlo = std::stof(split(tmp_str, "\"")[1]);

            // 0.000
            ylo = std::stof(tokens[10]);

            // 10.000"
            tmp_str = tokens[11];
            zlo = std::stof(split(tmp_str, "\"")[0]);

            has_boundaries = true;
            axis_lengths = {lx, ly, lz};

            continue;
        } else if (tokens.size() >= 4 && tokens.size() <= 5) {
            if (tokens.size() == 4) {
                atom_count++;
                atom_id = atom_count;
                given_type = tokens[0];
                x = std::stof(tokens[1]);
                y = std::stof(tokens[2]);
                z = std::stof(tokens[3]);
            } else if (tokens.size() == 5) {
                atom_count++;
                atom_id = std::stoi(tokens[0]);
                given_type = tokens[1];
                x = std::stof(tokens[2]);
                y = std::stof(tokens[3]);
                z = std::stof(tokens[4]);
            } else {
                throw std::runtime_error("Invalid atom line in xyz file: " + line);
            }

            // In xyz file, atom may have type of int (6...) or string (C...)
            if (!can_convert_to_int(given_type)) {
                // In case of string type, give a int type for atom.
                // Notice that must get the types correct!
                try {
                    if (type_stoi.find(given_type) == type_stoi.end()) {
                        type_int = type_stoi.size();
                        type_itos[type_int] = given_type;
                        type_stoi[given_type] = type_int;
                        itypes++;
                    } else {
                        type_int = type_stoi[given_type];
                    }
                    type_str = given_type;
                } catch (const std::invalid_argument& e) {
                    throw std::runtime_error("Invalid element type in xyzfile: " + given_type);
                }
            } else {
                // In case of int type
                try {
                    type_int = std::stoi(given_type);
                    type_str = type_itos[type_int];
                } catch (const std::invalid_argument& e) {
                    throw std::runtime_error("Invalid element type in xyzfile: " + given_type);
                }
            }

            if (has_boundaries) {
                x = x - xlo;
                y = y - ylo;
                z = z - zlo;
                x = fmod(x, lx);
                y = fmod(y, ly);
                z = fmod(z, lz);
                x = x < 0 ? x + lx : x;
                y = y < 0 ? y + ly : y;
                z = z < 0 ? z + lz : z;
            } else {
                x = x - xlo;
                y = y - ylo;
                z = z - zlo;
            }

            Atom* atom = new Atom(atom_id, type_int, {x, y, z}, type_str);
            atoms.push_back(atom);

            // When there's no atom numbers line in xyz file, iatoms = 0, and
            // first atom id = 1.
            if (atom_count == iatoms) {
                break;
            }
        }
    }
}

void System::load_lammpstrj(std::ifstream& file) {
    // load a frame from lammpstrj, create a system instance.
    // when finish all atoms in this frame, stop.
    std::string line;
    // For rectangular box, maybe support triclinc box later.
    // xlo, ylo, zlo, xhi, yhi, zhi
    std::vector<float> bounds(6, 0.0f);
    float xlo = 0.0f, xhi = 0.0f, ylo = 0.0f, yhi = 0.0f, zlo = 0.0f, zhi = 0.0f, lx = 0.0f, ly = 0.0f, lz = 0.0f;
    float x, y, z;

    // When lammpstrj style set to xs(,ys,zs), file uses relative coord.
    bool is_relative_coord = false;
    int atom_id_column = 0;
    int atom_type_column = 1;
    int atom_x_column = 2;
    int atom_y_column = 3;
    int atom_z_column = 4;

    bool read_timestep = false;
    bool read_natoms = false;
    bool read_box = false;
    int box_dim = 0;
    bool read_atoms = false;
    int atoms_count = 0;

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file");
    }

    while (getline(file, line)) {
        std::vector<std::string> tokens = split_by_space(line);
        if ((tokens.size() == 0) or (tokens[0] == "#")) {
            continue;
        } else if (tokens[0] == "ITEM:") {
            if (tokens[1] == "TIMESTEP") {
                read_timestep = true;
            }
            if (tokens[1] == "NUMBER") {
                read_natoms = true;
            }
            if (tokens[1] == "BOX") {
                read_box = true;
            }
            if (tokens[1] == "ATOMS") {
                if (std::find(tokens.begin(), tokens.end(), "xs") != tokens.end()) {
                    is_relative_coord = true;
                    atom_x_column = std::find(tokens.begin(), tokens.end(), "xs") - tokens.begin() - 2;
                    atom_y_column = std::find(tokens.begin(), tokens.end(), "ys") - tokens.begin() - 2;
                    atom_z_column = std::find(tokens.begin(), tokens.end(), "zs") - tokens.begin() - 2;
                } else {
                    // find index of id, type, x, y, z
                    atom_x_column = std::find(tokens.begin(), tokens.end(), "x") - tokens.begin() - 2;
                    atom_y_column = std::find(tokens.begin(), tokens.end(), "y") - tokens.begin() - 2;
                    atom_z_column = std::find(tokens.begin(), tokens.end(), "z") - tokens.begin() - 2;
                }
                atom_id_column = std::find(tokens.begin(), tokens.end(), "id") - tokens.begin() - 2;
                atom_type_column = std::find(tokens.begin(), tokens.end(), "type") - tokens.begin() - 2;

                static bool printed_column_info = false;
                if (!printed_column_info) {
                    fmt::print("Column index of lammpstrj file:id: {}, type: {}, x: {}, y: {}, z: {}\n", atom_id_column,
                               atom_type_column, atom_x_column, atom_y_column, atom_z_column);
                    printed_column_info = true;
                }

                read_atoms = true;
            }
        } else if (read_timestep) {
            itimestep = std::stoi(tokens[0]);
            read_timestep = false;
        } else if (read_natoms) {
            iatoms = std::stoi(tokens[0]);
            atoms.reserve(iatoms);
            bonds.reserve(iatoms * 3);
            molecules.reserve(iatoms / 2);
            read_natoms = false;
        } else if (read_box) {
            // For rectangular box, maybe support triclinc box later.
            if (box_dim == 0) {
                xlo = std::stof(tokens[0]);
                xhi = std::stof(tokens[1]);
                lx = xhi - xlo;
            } else if (box_dim == 1) {
                ylo = std::stof(tokens[0]);
                yhi = std::stof(tokens[1]);
                ly = yhi - ylo;
            } else if (box_dim == 2) {
                zlo = std::stof(tokens[0]);
                zhi = std::stof(tokens[1]);
                lz = zhi - zlo;
            }

            box_dim++;
            if (box_dim > 2) {
                bounds = {xlo, ylo, zlo, xhi, yhi, zhi};
                this->has_boundaries = true;
                axis_lengths = {lx, ly, lz};
                read_box = false;
            }
        } else if (read_atoms) {
            if (!has_boundaries) {
                throw std::runtime_error("lammstrj file must have boundaries!");
            }
            // Note: assume that the atom card style is "id type x y z" or "id
            // type xs ys zs" types must be set before
            int id = std::stoi(tokens[atom_id_column]);
            int type_i = std::stoi(tokens[atom_type_column]);
            std::string type_s = type_itos[type_i];

            // wrap and transform into {0, lx, 0, ly, 0, lz} box.

            if (!is_relative_coord) {
                x = std::stof(tokens[atom_x_column]) - xlo;
                y = std::stof(tokens[atom_y_column]) - ylo;
                z = std::stof(tokens[atom_z_column]) - zlo;
                x = fmod(x, lx);
                y = fmod(y, ly);
                z = fmod(z, lz);
                x = x < 0 ? x + lx : x;
                y = y < 0 ? y + ly : y;
                z = z < 0 ? z + lz : z;
            } else {
                x = std::stof(tokens[atom_x_column]) * lx;
                y = std::stof(tokens[atom_y_column]) * ly;
                z = std::stof(tokens[atom_z_column]) * lz;
            }

            Atom* atom = new Atom(id, type_i, {x, y, z}, type_s);
            atoms.push_back(atom);
            atoms_count++;

            // End reading.
            if (atoms_count == iatoms) {
                break;
            }
        }
    }
}