#include "fmt/format.h"
#include "string_tools.h"
#include "universe.h"

void System::load_xyz(std::ifstream& file) {
    std::string line;
    std::string delim = " ";
    int atom_id = 0;
    bool skip = false;

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file");
    }

    while (getline(file, line)) {
        std::vector<std::string> tokens = split_by_space(line);

        if (skip) {
            skip = false;
            continue;
        } else if ((tokens.size() == 0) or (tokens[0] == "#")) {
            continue;
        }
        // Atom numbers line
        else if (tokens.size() == 1 && can_convert_to_int(tokens[0])) {
            // System::iatoms
            iatoms = std::stoi(tokens[0]);
            atoms.reserve(iatoms);
            bonds.reserve(iatoms * 3);
            molecules.reserve(iatoms / 2);

            skip = true;
            continue;
        } else if (tokens.size() >= 4) {
            int type_i;
            atom_id++;
            std::string type_s;
            std::vector<float> coord;

            // In xyz file, atom may have type of int (6...) or string (C...)
            if (!can_convert_to_int(tokens[0])) {
                // In case of string type, give a int type for atom.
                type_s = tokens[0];
                if (type_stoi.contains(type_s)) {
                    type_i = type_stoi[type_s];
                } else {
                    type_i = ++itypes;
                    type_stoi[type_s] = type_i;
                    type_itos[type_i] = type_s;
                }
            } else {
                // In case of int type
                type_i = std::stoi(tokens[0]);
                if (!type_itos.contains(type_i)) {
                    type_s = fmt::format("E_{}_", type_i);
                    type_itos[type_i] = type_s;
                    type_stoi[type_s] = type_i;
                }
            }

            coord = {std::stof(tokens[1]), std::stof(tokens[2]), std::stof(tokens[3])};
            Atom* atom = new Atom(atom_id, type_i, coord, type_s);
            atoms.push_back(atom);

            // When there's no atom numbers line in xyz file, iatoms = 0, and
            // first atom id = 1.
            if (atom_id == iatoms) [[unlikely]] {
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
    float xlo, xhi, ylo, yhi, zlo, zhi, lx, ly, lz;

    // When lammpstrj style set to xs(,ys,zs), file uses relative coord.
    bool is_relative_coord = false;
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
                has_boundaries = true;
                axis_lengths = {lx, ly, lz};
                read_box = false;
            }
        } else if (read_atoms) {
            // Note: assume that the atom card style is "id type x y z" or "id
            // type xs ys zs" types must be set before
            int id = std::stoi(tokens[0]);
            int type_i = std::stoi(tokens[1]);
            std::string type_s = type_itos[type_i];
            std::vector<float> coord;
            float x, y, z;

            // wrap and transform into {0, lx, 0, ly, 0, lz} box.

            if (!is_relative_coord) {
                x = std::stof(tokens[2]);
                y = std::stof(tokens[3]);
                z = std::stof(tokens[4]);
            } else {
                x = std::stof(tokens[2]) * lx;
                y = std::stof(tokens[3]) * ly;
                z = std::stof(tokens[4]) * lz;
            }

            if (has_boundaries) {
                x = fmod(x, lx);
                y = fmod(y, ly);
                z = fmod(z, lz);
                x = x < 0 ? x + lx : x;
                y = y < 0 ? y + ly : y;
                z = z < 0 ? z + lz : z;
            }

            Atom* atom = new Atom(id, type_i, {x, y, z}, type_s);
            atoms.push_back(atom);
            atoms_count++;

            // End reading.
            if (atoms_count == iatoms) [[unlikely]] {
                break;
            }
        }
    }
}