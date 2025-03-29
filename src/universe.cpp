#include "defines.h"
#include "string_tools.h"

#include "universe.h"
#include "molecule.h"
#include "atom.h"
#include "reax_species.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_set>

#include "fmt/format.h"

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

void System::build_bonds_by_radius(const float& rvdw_scale)
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
				bond_radius[pair_ij] = 0.5f * (atomic_radius[type_i] + atomic_radius[type_j]) * rvdw_scale;
			}
			else if (type_itos[type_i] == "X" or type_itos[type_j] == "X")
			{
				bond_radius[pair_ij] = 0.0f;
			}
			else
			{
				bond_radius[pair_ij] = 0.5f * (atomic_radius[type_i] + atomic_radius[type_j]) * rvdw_scale;
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

void Universe::process_traj(const std::string& file_path, const std::vector<std::string>& type_names, const float& rvdw_scale)
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
		curr_sys->build_bonds_by_radius(rvdw_scale);
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
		// Update reaction flow.
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

