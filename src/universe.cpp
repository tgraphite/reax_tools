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

		// Species analysis. Only depend on molecule formulas.
		// TODO: Depend on topology.
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

