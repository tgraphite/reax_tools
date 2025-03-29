#pragma once
#include "defines.h"

#include <string>
#include <vector>
#include <map>
#include <fstream>

class ReaxSpecies
{
private:
	std::string file_path;
	std::string save_path;
	std::ifstream file;
	std::ofstream save;
	std::vector<std::string> order;
	std::vector<std::string> formulas;
	std::map<std::string, std::vector<float>> formulas_nums;
	std::vector<std::vector<std::string>> all_frame_formulas;
	int nframes = 0;

public:
	// When read from Lammps reax/c/bonds file.
	ReaxSpecies(const std::string& file_path, const std::string& save_path,
		const std::vector<std::string>& order = sorted_elements);
	ReaxSpecies(const std::string& file_path, const std::vector<std::string>& order = sorted_elements);
	ReaxSpecies(const std::vector<std::string>& order = sorted_elements);

	~ReaxSpecies();

	void get_formulas();
	void get_nums();
	void import_frame_formulas(const std::vector<std::string>& frame_formulas);
	void analyze_frame_formulas();

	void rename_all_formulas();
	void merge_formulas(const std::vector<std::string>& formulas, const std::string& new_formula);
	void scale_formula(const std::string& formula, const float& k);
	void rescale_all_by_element(const std::string& target_element);
	void merge_by_element(const std::string& target_element, const std::vector<int>& ranges, bool recalc = false);

	void show_nums();
	void brief_report();

	// TODO
	void save_file();
	void save_file(const std::string& raw_file_path);
};
