#include "defines.h"
#include "string_tools.h"

#include "universe.h"
#include "molecule.h"
// #include "mol_topo.h"
#include "atom.h"
#include "reax_species.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_set>

#include "fmt/format.h"

class ReaxFlow;  // 前向声明

void System::search_neigh(const float& radius, const int& max_neigh)
{
	KD_tree kd_tree;

	for (auto& atom : atoms)
	{
		kd_tree.insert(atom);
	}

	for (auto& atom_data : atoms)
	{
		if (atom_data->neighs.size() < max_neigh)
		{
			std::vector<std::shared_ptr<Atom>> neighbors;
			kd_tree.find_neighbors(atom_data, radius, neighbors);
			for (std::shared_ptr<Atom> neighbor : neighbors)
			{
				atom_data->neighs.push_back(neighbor);
				if (atom_data->neighs.size() >= max_neigh)
				{
					break;
				}
			}
		}
	}
}

void System::build_bonds_by_radius(const float& tolerance)
{
	std::map<int, float> atomic_radius;
	/// When using default atomic radius in constant.h
	for (const auto& pair : type_stoi)
	{
		std::string typ_s = pair.first;
		int typ_i = pair.second;

		// pair.first: string type, pair.second: int type
		if (default_atomic_radius.contains(typ_s))
			atomic_radius[typ_i] = default_atomic_radius.at(typ_s);
		else
			atomic_radius[typ_i] = default_atomic_radius.at("X");
	}

	for (int type_i = 1; type_i <= itypes; type_i++)
	{
		for (int type_j = 1; type_j <= itypes; type_j++)
		{
			std::pair<int, int> pair_ij = { type_i, type_j };
			if (type_itos[type_i] == "H" && type_itos[type_j] == "H")
			{
				bond_radius[pair_ij] = 0.5f * (atomic_radius[type_i] + atomic_radius[type_j]) * tolerance;
			}
			else if (type_itos[type_i] == "X" or type_itos[type_j] == "X")
			{
				bond_radius[pair_ij] = 0.0f;
			}
			else
			{
				bond_radius[pair_ij] = 0.5f * (atomic_radius[type_i] + atomic_radius[type_j]) * tolerance;
			}
		}
	}

	float bond_r;
	float bond_sq;
	float dist_sq;
	float dx, dy, dz;

	for (auto& atom : atoms)
	{
		for (auto& neigh : atom->neighs)
		{
			if (neigh == atom)
			{
				continue;
			}
			if (std::find(neigh->bonded_atoms.begin(), neigh->bonded_atoms.end(), atom) != neigh->bonded_atoms.end())
			{
				continue;
			}

			dx = atom->coord[0] - neigh->coord[0];
			dy = atom->coord[1] - neigh->coord[1];
			dz = atom->coord[2] - neigh->coord[2];

			dist_sq = dx * dx + dy * dy + dz * dz;

			std::pair<int, int> id_ij = { atom->type_id, neigh->type_id };
			bond_r = bond_radius[id_ij];
			bond_sq = bond_r * bond_r;

			if (dist_sq < bond_sq)
			{
				std::pair<std::shared_ptr<Atom>, std::shared_ptr<Atom>> pair_ij = { atom, neigh };
				auto bond = std::make_shared<Bond>(atom, neigh);

				bonds.push_back(bond);
				atom->bonds.push_back(bond);
				neigh->bonds.push_back(bond);
				atom->bonded_atoms.push_back(neigh);
				neigh->bonded_atoms.push_back(atom);
			}
		}
	}
}

void System::build_molecules()
{
	std::set<std::shared_ptr<Atom>> visited;
	int new_mol_id = 0;

	for (auto& atom : atoms)
	{
		if (visited.find(atom) == visited.end()) // If not visited
		{
			new_mol_id++;
			auto new_mol = std::make_shared<Molecule>(new_mol_id);
			dfs(atom, visited, new_mol);
			molecules.push_back(new_mol);
		}
	}

	for (auto& molecule : molecules)
	{
		molecule->update_formula();
		molecule->update_topo();
		// uint32_t hash = moltopo_bank->include_mol(molecule);
		// molecule->hash = hash;
	}
}

void System::dfs(std::shared_ptr<Atom>& atom, std::set<std::shared_ptr<Atom>>& visited, std::shared_ptr<Molecule>& cur_mol)
{
	visited.insert(atom);
	cur_mol->insert(atom);

	for (auto& bond : atom->bonds)
	{
		cur_mol->insert(bond);
		// One of atom_i and atom_j is the atom itself.
		if (visited.find(bond->atom_i) == visited.end())
		{
			dfs(bond->atom_i, visited, cur_mol);
			continue;
		}
		if (visited.find(bond->atom_j) == visited.end())
		{
			dfs(bond->atom_j, visited, cur_mol);
		}
	}
}



/*KD_tree*/
void KD_tree::insert(std::shared_ptr<Atom>& atom_data)
{
	root = insert_rec(root, atom_data, 0);
}

KD_tree::node* KD_tree::insert_rec(KD_tree::node* node_ptr, std::shared_ptr<Atom>& atom_data, const int& depth)
{
	if (!node_ptr)
		return new node(atom_data, depth % 3);

	int axis = node_ptr->axis;
	if (atom_data->coord[axis] < node_ptr->atom_data->coord[axis])
		node_ptr->left = insert_rec(node_ptr->left, atom_data, depth + 1);
	else
		node_ptr->right = insert_rec(node_ptr->right, atom_data, depth + 1);

	return node_ptr;
}

void KD_tree::find_neighbors(const std::shared_ptr<Atom>& atom_data, const float& radius, std::vector<std::shared_ptr<Atom>>& neighbors)
{
	find_neighbors_rec(root, atom_data, radius, neighbors);
}

void KD_tree::find_neighbors_rec(KD_tree::node* node_ptr, const std::shared_ptr<Atom>& atom_data, const float& radius, std::vector<std::shared_ptr<Atom>>& neighbors)
{
	if (!node_ptr)
		return;

	float dx = node_ptr->atom_data->coord[0] - atom_data->coord[0];
	float dy = node_ptr->atom_data->coord[1] - atom_data->coord[1];
	float dz = node_ptr->atom_data->coord[2] - atom_data->coord[2];
	float dist_sq = dx * dx + dy * dy + dz * dz;

	if (dist_sq < radius * radius)
		neighbors.push_back(node_ptr->atom_data);

	int axis = node_ptr->axis;
	float delta = atom_data->coord[axis] - node_ptr->atom_data->coord[axis];
	float delta_sq = delta * delta;
	float radius_sq = radius * radius;
	if (delta < 0)
	{
		find_neighbors_rec(node_ptr->left, atom_data, radius, neighbors);
		if (delta_sq < radius_sq)
			find_neighbors_rec(node_ptr->right, atom_data, radius, neighbors);
	}
	else
	{
		find_neighbors_rec(node_ptr->right, atom_data, radius, neighbors);
		if (delta_sq < radius_sq)
			find_neighbors_rec(node_ptr->left, atom_data, radius, neighbors);
	}
}

void KD_tree::clear(KD_tree::node* node_ptr)
{
	if (node_ptr)
	{
		clear(node_ptr->left);
		clear(node_ptr->right);
		delete node_ptr;
	}
}

KD_tree::~KD_tree()
{
	clear(root);
}

System::System() {}

System::~System()
{
	// Have to clear all the low->high pointers first, otherwise will cause cycle reference and memory leak.
	for (auto& atom : atoms)
	{
		atom->bonds.clear();
		atom->neighs.clear();
		atom->bonded_atoms.clear();
	}
	atoms.clear();
	molecules.clear();
	dihedrals.clear();
	angles.clear();
	bonds.clear();
}

void System::set_types(const std::vector<std::string>& type_names)
{
	itypes = type_names.size();
	int type_id = 0;
	for (auto& name : type_names)
	{
		type_id++;
		type_itos[type_id] = name;
		type_stoi[name] = type_id;
	}
}

// void System::link_moltopo_bank(const std::shared_ptr<MolTopoBank>& _moltopo_bank)
// {
// 	moltopo_bank = _moltopo_bank;
// }

void System::load_xyz(std::ifstream& file)
{
	std::string line;
	std::string delim = " ";
	int atom_id = 0;
	bool skip = false;

	if (!file.is_open())
	{
		throw std::runtime_error("Failed to open file");
	}

	while (getline(file, line))
	{
		std::vector<std::string> tokens = split_by_space(line);

		if (skip)
		{
			skip = false;
			continue;
		}
		else if ((tokens.size() == 0) or (tokens[0] == "#"))
		{
			continue;
		}
		// Atom numbers line
		else if (tokens.size() == 1 && can_convert_to_int(tokens[0]))
		{
			// System::iatoms
			iatoms = std::stoi(tokens[0]);
			atoms.reserve(iatoms);
			bonds.reserve(iatoms * 3);
			angles.reserve(iatoms * 5);
			dihedrals.reserve(iatoms * 8);
			molecules.reserve(iatoms / 2);

			skip = true;
			continue;
		}
		else if (tokens.size() >= 4)
		{
			int type_i;
			atom_id++;
			std::string type_s;
			std::vector<float> coord;

			// In xyz file, atom may have type of int (1,2...) or string (C, H...)
			if (!can_convert_to_int(tokens[0]))
			{
				// In case of string type, give a int type for atom.
				type_s = tokens[0];
				if (type_stoi.contains(type_s))
				{
					type_i = type_stoi[type_s];
				}
				else
				{
					type_i = ++itypes;
					type_stoi[type_s] = type_i;
					type_itos[type_i] = type_s;
				}
			}
			else
			{
				// In case of int type
				type_i = std::stoi(tokens[0]);
				if (!type_itos.contains(type_i))
				{
					type_s = fmt::format("E_{}_", type_i);
					type_itos[type_i] = type_s;
					type_stoi[type_s] = type_i;
				}
			}

			coord = { std::stof(tokens[1]), std::stof(tokens[2]), std::stof(tokens[3]) };
			std::shared_ptr<Atom> atom = std::make_shared<Atom>(atom_id, type_i, coord, type_s);
			atoms.push_back(atom);

			// When there's no atom numbers line in xyz file, iatoms = 0, and first atom id = 1.
			if (atom_id == iatoms) [[unlikely]]
			{
				break;
			}
		}
	}
}

void System::load_lammpstrj(std::ifstream& file)
{
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

	if (!file.is_open())
	{
		throw std::runtime_error("Failed to open file");
	}

	while (getline(file, line))
	{
		std::vector<std::string> tokens = split_by_space(line);
		if ((tokens.size() == 0) or (tokens[0] == "#"))
		{
			continue;
		}
		else if (tokens[0] == "ITEM:")
		{
			if (tokens[1] == "TIMESTEP")
			{
				read_timestep = true;
			}
			if (tokens[1] == "NUMBER")
			{
				read_natoms = true;
			}
			if (tokens[1] == "BOX")
			{
				read_box = true;
			}
			if (tokens[1] == "ATOMS")
			{
				if (std::find(tokens.begin(), tokens.end(), "xs") != tokens.end())
				{
					is_relative_coord = true;
				}
				read_atoms = true;
			}
		}
		else if (read_timestep)
		{
			itimestep = std::stoi(tokens[0]);
			read_timestep = false;
		}
		else if (read_natoms)
		{
			iatoms = std::stoi(tokens[0]);
			atoms.reserve(iatoms);
			bonds.reserve(iatoms * 3);
			angles.reserve(iatoms * 5);
			dihedrals.reserve(iatoms * 8);
			molecules.reserve(iatoms / 2);
			read_natoms = false;
		}
		else if (read_box)
		{
			// For rectangular box, maybe support triclinc box later.
			if (box_dim == 0)
			{
				xlo = std::stof(tokens[0]);
				xhi = std::stof(tokens[1]);
				lx = xhi - xlo;
			}
			else if (box_dim == 1)
			{
				ylo = std::stof(tokens[0]);
				yhi = std::stof(tokens[1]);
				ly = yhi - ylo;
			}
			else if (box_dim == 2)
			{
				zlo = std::stof(tokens[0]);
				zhi = std::stof(tokens[1]);
				lz = zhi - zlo;
			}

			box_dim++;
			if (box_dim > 2)
			{
				bounds = { xlo, ylo, zlo, xhi, yhi, zhi };
				read_box = false;
			}
		}
		else if (read_atoms)
		{
			// Note: assume that the atom card style is "id type x y z" or "id type xs ys zs"
			// types must be set before
			int id = std::stoi(tokens[0]);
			int type_i = std::stoi(tokens[1]);
			std::string type_s = type_itos[type_i];
			std::vector<float> coord;

			if (!is_relative_coord)
			{
				coord = { std::stof(tokens[2]), std::stof(tokens[3]), std::stof(tokens[4]) };
			}
			else
			{
				coord =
				{
					xlo + std::stof(tokens[2]) * lx,
					ylo + std::stof(tokens[3]) * ly,
					zlo + std::stof(tokens[4]) * lz,
				};
			}

			std::shared_ptr<Atom> atom = std::make_shared<Atom>(id, type_i, coord, type_s);
			atoms.push_back(atom);
			atoms_count++;

			// End reading.
			if (atoms_count == iatoms) [[unlikely]]
			{
				break;
			}
		}
	}
}

void System::summary()
{
	fprintf(stdout, "Atoms: %.6d Bonds: %.6d Mols: %.6d\n", atoms.size(), bonds.size(), molecules.size());
	std::cout << "Atoms: " << atoms.size() << std::endl;
	print_info(atoms, 5, true);
	std::cout << "Bonds: " << bonds.size() << std::endl;
	print_info(bonds, 5, true);
	std::cout << "Mols: " << molecules.size() << std::endl;
	print_info(molecules, 5, true);
}

void System::basic_info()
{
	std::cout << "Atoms: " << atoms.size() << " Bonds: " << bonds.size() << " Mols: " << molecules.size() << std::endl;
}

Universe::Universe() {};

Universe::~Universe() {
	prev_sys.reset();
	curr_sys.reset();
}

void Universe::flush()
{
	prev_sys = curr_sys;
	curr_sys.reset();
}

// Tick-tock reading two frames and compare.
// args:
// file_path: string, path-like;
// type_names: vector<string> element symbols {"C", "H", ...}
// reax_speices: empty ReaxSpecies instance for get data from coordinate analysis.
void Universe::process_traj(const std::string& file_path, const std::vector<std::string>& type_names, const float& tolerance)
{
	std::ifstream file(file_path);
	int curr_frame = 1;

	while (file.is_open() and !file.eof())
	{
		flush();
		curr_sys = std::make_shared<System>();
		curr_sys->set_types(type_names);
		// curr_sys->link_moltopo_bank(moltopo_bank);

		if (file_path.ends_with(".lammpstrj"))
			curr_sys->load_lammpstrj(file);
		else if (file_path.ends_with(".xyz"))
			curr_sys->load_xyz(file);

		curr_sys->search_neigh(3.0, 10);
		curr_sys->build_bonds_by_radius(tolerance);
		curr_sys->build_molecules();

		if (curr_sys->atoms.size() == 0)
			continue;

		// Species analysis.
		std::vector<std::string> frame_formulas;
		frame_formulas.reserve(curr_sys->molecules.size());
		for (auto& mol : curr_sys->molecules)
		{
			frame_formulas.push_back(mol->formula);
		}
		reax_species->import_frame_formulas(frame_formulas);

		// curr_sys->summary();
		// Update evolution.
		if (prev_sys && curr_sys) {
			update_reax_flow(prev_sys, curr_sys, curr_frame);
		}
		// compare prev_sys and curr_sys
		std::cout << "Frame: " << curr_frame << " ";
		curr_frame++;
		curr_sys->basic_info();
	}
	reax_species->analyze_frame_formulas();
}

void Universe::update_reax_flow(const std::shared_ptr<System>& prev_sys, 
                               const std::shared_ptr<System>& curr_sys,
                               const int& curr_frame) {
								
    for (auto& prev_mol : prev_sys->molecules) {
        // Ignore single atom molecule.
        if (prev_mol->atom_ids.size() == 1)
            continue;
            
        // Find the most similar molecule in current frame.
		
        std::shared_ptr<Molecule> best_match = nullptr;
        float best_similarity = 0.0f;
        
        std::vector<int> intersection;
        std::vector<int> union_set;
        
        // Optimization: stop searching once a similarity > 0.5 is found.
        bool found_match = false;
        
        for (auto& curr_mol : curr_sys->molecules) {
            // Ignore single atom molecule.
            if (curr_mol->atom_ids.size() == 1)
                continue;
                
            // If the formula is the same, consider it the same molecule, skip.
            if (prev_mol->formula == curr_mol->formula)
                continue;  
                
            // Quick filter: skip if the size difference is too large.
            if ((curr_mol->atom_ids.size() / prev_mol->atom_ids.size() >= 2) ||
                (prev_mol->atom_ids.size() / curr_mol->atom_ids.size() >= 2))
                continue;
                
            // Calculate similarity: intersection / union.
            intersection.clear();
            union_set.clear();
            
            std::set_intersection(prev_mol->atom_ids.begin(), prev_mol->atom_ids.end(),
                curr_mol->atom_ids.begin(), curr_mol->atom_ids.end(), 
                back_inserter(intersection));
            
            // Quick check: if intersection is empty, similarity is 0.
            if (intersection.empty())
                continue;
                
            // Quick check: if intersection size equals prev_mol size and curr_mol size, then they are the same molecule.
            if (intersection.size() == prev_mol->atom_ids.size() && 
                intersection.size() == curr_mol->atom_ids.size())
                continue;  // Same molecule, not reaction.
                
            std::set_union(prev_mol->atom_ids.begin(), prev_mol->atom_ids.end(),
                curr_mol->atom_ids.begin(), curr_mol->atom_ids.end(), 
                back_inserter(union_set));
                
            float similarity = float(intersection.size()) / float(union_set.size());
            
            // If similarity exceeds the threshold, record and stop searching.
            if (similarity >= 0.5) {
                best_match = curr_mol;
                best_similarity = similarity;
                found_match = true;
                break;  // If a good match is found, stop searching.
            }
            
            // If similarity exceeds the threshold, record and stop searching.
            else if (similarity > best_similarity) {
                best_match = curr_mol;
                best_similarity = similarity;
            }
        }
        
        // If found a match and not the same molecule (similarity between 0.5 and 1.0), record the reaction.
        if (best_match && best_similarity >= 0.5 && best_similarity < 1.0) {
            reax_flow->add_reaction(curr_frame, prev_mol, best_match);
        }
    }
}