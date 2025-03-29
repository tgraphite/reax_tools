#include "universe.h"
#include "atom.h"
#include "molecule.h"
#include "fmt/format.h"
#include "defines.h"

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
    fmt::print("Atoms: {}, Bonds: {}, Mols: {}\n", atoms.size(), bonds.size(), molecules.size());
}

void System::basic_info()
{
	fmt::print("Atoms: {}, Bonds: {}, Mols: {}\n", atoms.size(), bonds.size(), molecules.size());
}

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