#include "reax_species.h"

#include <iostream>
#include <mutex>
#include <string>
#include <unordered_set>
#include <vector>

#include "fmt/format.h"
#include "string_tools.h"
#include "vec_algorithms.h"

// For read data from Lammps species file.
ReaxSpecies::ReaxSpecies(const std::string &file_path) : file_path(file_path) {
    file.open(file_path);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << file_path << std::endl << "Exit." << std::endl;
        exit(EXIT_FAILURE);
    }

    get_formulas();
    get_nums();
}

// For import data from class Universe.
ReaxSpecies::ReaxSpecies() {}

ReaxSpecies::~ReaxSpecies() {
    if (file.is_open()) file.close();
}

void ReaxSpecies::get_formulas() {
    std::string line;
    std::unordered_set<std::string> formulas_set;

    while (getline(file, line)) {
        if (line[0] == '#') {
            nframes++;
            std::vector<std::string> parts = split_by_space(line);

            for (std::string formula : parts) {
                if (formula == "#" || formula == "Timestep" || formula == "No_Moles" || formula == "No_Specs") {
                    continue;
                } else {
                    formulas_set.insert(formula);
                }
            }
        }
    }
    file.clear();
    file.seekg(0);

    formulas = std::vector<std::string>(formulas_set.begin(), formulas_set.end());
    for (const std::string &formula : formulas) {
        formulas_nums[formula] = std::vector<float>(nframes, 0.0);
    }
    return;
}

void ReaxSpecies::get_nums() {
    std::string line;
    std::vector<std::string> cur_formulas;
    int cur_frame = 0;

    while (getline(file, line)) {
        int cur_step, cur_mols, cur_spec;
        std::vector<std::string> parts = split_by_space(line);

        if (parts[0] == "#") {
            cur_formulas.clear();

            for (size_t i = 4; i < parts.size(); i++) {
                // std::string formula = rename_species(parts[i]);
                std::string formula = parts[i];
                cur_formulas.push_back(formula);
            }
        } else {
            cur_step = std::stoi(parts[0]);
            cur_mols = std::stoi(parts[1]);
            cur_spec = std::stoi(parts[2]);

            for (size_t i = 0; i < cur_formulas.size(); i++) {
                std::string cur_formula = cur_formulas[i];
                float cur_num = std::stof(parts[i + 3]);
                formulas_nums[cur_formula][cur_frame] = cur_num;
            }
            cur_frame++;
        }
    }
}

void ReaxSpecies::rename_all_formulas(const std::vector<std::string> &order) {
    // Create a temporary map to store formulas that need to be renamed.
    std::map<std::string, std::vector<float>> to_insert;
    std::vector<std::string> to_erase;

    // First step: collect all formulas that need to be renamed.
    for (const auto &pair : formulas_nums) {
        std::string old_formula = pair.first;
        std::vector<float> nums = pair.second;

        if (starts_with(old_formula, "grp_")) {
            continue;
        }

        std::string new_formula = rename_formula(old_formula, order);

        if (!(new_formula == old_formula)) {
            to_insert[new_formula] = nums;
            to_erase.push_back(old_formula);
        }
    }

    // Second step: delete old formulas and add new formulas.
    for (const auto &formula : to_erase) {
        auto it = formulas_nums.find(formula);
        if (it != formulas_nums.end()) {
            formulas_nums.erase(it);
        }
    }

    // Third step: add new formulas.
    for (const auto &pair : to_insert) {
        formulas_nums[pair.first] = pair.second;
    }
}

// User accessable.
void ReaxSpecies::merge_formulas(const std::vector<std::string> &formulas, const std::string &new_formula) {
    std::vector<float> new_nums(nframes, 0.0);
    for (const std::string &formula : formulas) {
        std::vector<float> num = formulas_nums[formula];
        new_nums = add_vectors(new_nums, num);
    }

    for (const std::string &formula : formulas) {
        formulas_nums.erase(formula);
    }

    formulas_nums[new_formula] = new_nums;
    return;
}

// User accessable.
void ReaxSpecies::scale_formula(const std::string &formula, const float &k) {
    std::vector<float> new_nums(nframes, 0.0);
    new_nums = vector_scale(formulas_nums[formula], k);
    formulas_nums[formula] = new_nums;
    return;
}

// User accessable.
void ReaxSpecies::rescale_all_by_element(const std::string &target_element) {
    for (auto &pair : formulas_nums) {
        std::string formula = pair.first;
        std::vector<float> &old_nums = pair.second;
        std::map<std::string, int> elements_weights = parse_formula(formula);

        if (starts_with(formula, "group_")) {
            continue;
        }

        for (auto &elem_weight : elements_weights) {
            std::string elem = elem_weight.first;
            float weight = float(elem_weight.second);

            if (target_element == elem) {
                std::vector<float> new_nums = vector_scale(old_nums, weight);
                formulas_nums[formula] = new_nums;
            }
        }
    }
    return;
}

// User accessable.
void ReaxSpecies::merge_by_element(const std::string &target_element, const std::vector<int> &ranges, bool rescale) {
    for (size_t i = 0; i < ranges.size(); i++) {
        int start;
        int end;
        std::string new_formula;

        if (i < ranges.size() - 1) {
            start = ranges[i];
            end = ranges[i + 1] - 1;
            new_formula = "grp_" + target_element + std::to_string(start) + "-" + std::to_string(end);
        } else {
            start = ranges[i];
            end = 10000;
            new_formula = "grp_>" + target_element + std::to_string(start);
        }

        std::vector<std::string> formulas_to_erase;
        std::vector<float> new_nums(nframes, 0.0);

        for (auto &pair : formulas_nums) {
            std::string old_formula = pair.first;
            std::vector<float> &old_nums = pair.second;

            if (starts_with(old_formula, "grp_")) {
                continue;
            }

            std::map<std::string, int> elements_weights = parse_formula(old_formula);

            for (auto &elem_weight : elements_weights) {
                std::string elem = elem_weight.first;
                float weight = float(elem_weight.second);

                if (target_element == elem && weight >= start && weight <= end) {
                    if (rescale) {
                        std::vector<float> weighed_nums = vector_scale(old_nums, weight);
                        new_nums = add_vectors(new_nums, weighed_nums);
                        formulas_nums[old_formula] = std::vector<float>(nframes, 0.0);
                        formulas_to_erase.push_back(old_formula);
                    } else {
                        new_nums = add_vectors(new_nums, old_nums);
                        formulas_nums[old_formula] = std::vector<float>(nframes, 0.0);
                        formulas_to_erase.push_back(old_formula);
                    }
                }
            }
        }

        formulas_nums[new_formula] = new_nums;
        std::cout << "Merge " << formulas_to_erase.size() << " formulas into " << new_formula << std::endl;

        for (const auto &formula : formulas_to_erase) {
            formulas_nums.erase(formula);
        }
    }
    return;
}

void ReaxSpecies::brief_report() {
    if (nframes < 40) {
        this->show_nums();
        return;
    }

    std::string header = fmt::format("{:<20s}{:>8s}{:>8s}{:>8s}{:>8s}", "formula", "begin", "mid", "end", "average");
    std::cout << header << std::endl;

    for (const auto &pair : formulas_nums) {
        std::string formula = pair.first;
        std::vector<float> nums = pair.second;
        size_t range = nframes / 10;
        size_t begin_l = 0;
        size_t mid_l = nframes * 5 / 10;
        size_t end_l = nframes * 9 / 10;

        float f_range = float(range);
        float average = 0.0f;
        float begin = 0.0f;
        float mid = 0.0f;
        float end = 0.0f;

        for (size_t i = 0; i < nframes; i++) {
            average += nums[i];
            if ((i >= begin_l) and (i < begin_l + range))
                begin += nums[i];
            else if ((i >= mid_l) and (i < mid_l + range))
                mid += nums[i];
            else if ((i >= end_l) and (i < end_l + range))
                end += nums[i];
        }

        average /= float(nframes);
        begin /= f_range;
        mid /= f_range;
        end /= f_range;

        std::string info = fmt::format("{:<20s}{:>8.2f}{:>8.2f}{:>8.2f}{:>8.2f}", formula, begin, mid, end, average);
        std::cout << info << std::endl;
    }
}

void ReaxSpecies::show_nums() {
    for (const auto &pair : formulas_nums) {
        std::cout << pair.first << " : ";
        if (nframes > 20) {
            for (size_t i = 0; i < 10; i++) {
                std::cout << pair.second[i] << " ";
            }
            std::cout << "...";
            for (size_t i = nframes - 10; i < nframes; i++) {
                std::cout << pair.second[i] << " ";
            }
        } else {
            for (size_t i = 0; i < nframes; i++) {
                std::cout << pair.second[i] << " ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    return;
}

// Get current frame formulas (std::map<std::string, int>) from class Universe.
void ReaxSpecies::import_frame_formulas(const std::vector<std::string> &formulas) {
    // 使用互斥锁保护共享资源
    std::lock_guard<std::mutex> lock(mutex_);

    all_frame_formulas.push_back(formulas);
}

// Analyze imported frame formulas and get standard formulas_nums after
// ReaxSpeices::import_frame_formulas all done.
void ReaxSpecies::analyze_frame_formulas() {
    // 使用互斥锁保护共享资源
    std::lock_guard<std::mutex> lock(mutex_);

    // all_formulas_set : any occured formula, non-repeat.
    std::unordered_set<std::string> all_formulas_set;
    nframes = all_frame_formulas.size();

    for (auto &frame_formulas_set : all_frame_formulas) {
        for (auto &formula : frame_formulas_set) {
            // If this formula not in all_formulas_set, append, otherwise skip.
            all_formulas_set.insert(formula);
        }
    }

    // Initialize formulas_nums for all formulas appeared in class Universe.
    // formulas_nums: map<string, vector<float>>
    formulas_nums.clear();
    for (const std::string &formula : all_formulas_set) {
        formulas_nums[formula] = std::vector<float>(nframes, 0.0);
    }

    for (size_t i = 0; i < nframes; i++) {
        for (auto &formula : all_frame_formulas[i]) {
            formulas_nums[formula][i] += 1.0f;
        }
    }
}

void ReaxSpecies::save_file() {
    save_path = file_path.substr(0, file_path.find_last_of(".")) + ".csv";

    FILE *file = fopen(save_path.c_str(), "w");
    if (file == NULL) {
        fprintf(stderr, "Cannot open file %s for writing.\n", save_path.c_str());
        return;
    }

    bool is_first = true;
    for (const auto &pair : formulas_nums) {
        if (is_first) {
            fprintf(file, "%s", pair.first.c_str());
            is_first = false;
        } else {
            fprintf(file, ",%s", pair.first.c_str());
        }
    }
    fprintf(file, "\n");

    for (size_t i = 0; i < nframes; i++) {
        is_first = true;
        for (const auto &pair : formulas_nums) {
            if (is_first) {
                fprintf(file, "%.0f", (pair.second)[i]);
                is_first = false;
            } else {
                fprintf(file, ",%.0f", (pair.second)[i]);
            }
        }
        fprintf(file, "\n");
    }

    fclose(file);
    fprintf(stdout, "Save file %s successfully.\n", save_path.c_str());
}

// When calling from Universe, use this function to set file_path and save.
// In that situation, the constructor did not get the raw file path, so we need
// to set it here.
void ReaxSpecies::save_file(const std::string &raw_file_path) {
    this->file_path = raw_file_path;
    save_file();
}