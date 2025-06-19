#include "reax_species.h"

#include <algorithm>
#include <iostream>
#include <numeric>
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

    // First step: collect all formulas that need to be renamed.

    std::string new_formula;
    std::string old_formula;
    std::vector<float> nums(nframes);
    std::map<std::string, std::vector<float>> new_formulas_nums;

    std::map<std::string, std::string> rename_map;

    for (const auto &pair : formulas_nums) {
        old_formula = pair.first;
        nums = pair.second;

        new_formula = rename_formula(old_formula, order);

        new_formulas_nums[new_formula] = nums;
    }

    formulas_nums.clear();
    formulas_nums = std::move(new_formulas_nums);
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
            new_formula = fmt::format("grp_{}{}-{}", target_element, start, end);
        } else {
            start = ranges[i];
            end = 10000;
            new_formula = fmt::format("grp_{}{}-max", target_element, start);
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
        fmt::print("Merge {} formulas into {}\n", formulas_to_erase.size(), new_formula);

        for (const auto &formula : formulas_to_erase) {
            formulas_nums.erase(formula);
        }
    }
    fmt::print("\n");
    return;
}

void ReaxSpecies::brief_report() {
    // if (nframes < 40) {
    //     this->show_nums();
    //     return;
    // }

    // sort formulas_nums by average in descending order
    std::vector<std::pair<std::string, float>> sorted_formulas_averages;
    float average = 0;
    for (const auto &pair : formulas_nums) {
        average = std::accumulate(pair.second.begin(), pair.second.end(), 0.0f) / nframes;
        sorted_formulas_averages.push_back(std::make_pair(pair.first, average));
    }
    std::sort(sorted_formulas_averages.begin(), sorted_formulas_averages.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });

    // get formulas going to report
    int max_species_print = 20;
    int min_frames_print = 10;
    if (sorted_formulas_averages.size() > max_species_print) {
        fmt::print("Note: To avoid too much screen output, only the top {} species are printed.\n", max_species_print);
        sorted_formulas_averages.resize(max_species_print);
    }

    if (nframes < min_frames_print) {
        fmt::print("Note: The number of frames is less than {}.\n", min_frames_print);
        fmt::print("{}\n", "=== Species Report ===");
        std::string header = fmt::format("{:<20s}{:>8s}{:>8s}{:>8s}", "formula", "begin", "mid", "end");
        fmt::print("{}\n", header);
        for (const auto &pair : sorted_formulas_averages) {
            std::string formula = pair.first;
            std::vector<float> nums = formulas_nums[formula];
            int mid_frame = nframes / 2;
            fmt::print("{:<20s}{:>8.2f}{:>8.2f}{:>8.2f}\n", formula, nums[0], nums[mid_frame], nums[nframes - 1]);
        }
    } else {
        fmt::print("{}\n", "=== Species Report ===");
        std::string header =
            fmt::format("{:<20s}{:>8s}{:>8s}{:>8s}{:>8s}", "formula", "begin", "mid", "end", "average");
        fmt::print("{}\n", header);
        for (const auto &pair : sorted_formulas_averages) {
            std::string formula = pair.first;
            std::vector<float> nums = formulas_nums[formula];
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
                if ((i >= begin_l) && (i < begin_l + range))
                    begin += nums[i];
                else if ((i >= mid_l) && (i < mid_l + range))
                    mid += nums[i];
                else if ((i >= end_l) && (i < end_l + range))
                    end += nums[i];
            }

            average /= float(nframes);
            begin /= f_range;
            mid /= f_range;
            end /= f_range;

            fmt::print("{:<20s}{:>8.2f}{:>8.2f}{:>8.2f}{:>8.2f}\n", formula, begin, mid, end, average);
        }
    }
    fmt::print("\n");
}

// Get current frame formulas (std::map<std::string, int>) from class Universe.
void ReaxSpecies::import_frame_formulas(int &frame_id, const std::vector<std::string> &formulas) {
    if (all_frame_formulas.size() < frame_id) {
        all_frame_formulas.resize(frame_id);
    }

    all_frame_formulas[frame_id - 1] = formulas;
}

// Analyze imported frame formulas and get standard formulas_nums after
// ReaxSpeices::import_frame_formulas all done.
void ReaxSpecies::analyze_frame_formulas() {
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

// General function, save file to given path.
void ReaxSpecies::save_file(const std::string &save_path) {
    FILE *file = fopen(save_path.c_str(), "w");
    if (file == NULL) {
        fmt::print(stderr, "Cannot open file {} for writing.\n", save_path);
        return;
    }

    bool is_first = true;
    for (const auto &pair : formulas_nums) {
        if (is_first) {
            fmt::print(file, "{}", pair.first);
            is_first = false;
        } else {
            fmt::print(file, ",{}", pair.first);
        }
    }
    fmt::print(file, "\n");

    for (size_t i = 0; i < nframes; i++) {
        is_first = true;
        for (const auto &pair : formulas_nums) {
            if (is_first) {
                fmt::print(file, "{}", (pair.second)[i]);
                is_first = false;
            } else {
                fmt::print(file, ",{}", (pair.second)[i]);
            }
        }
        fmt::print(file, "\n");
    }

    fclose(file);
    fmt::print(stdout, "Save species file {} successfully.\n\n", save_path);
}

// For traj mode, call from universe.
void ReaxSpecies::save_file_to_dir(const std::string &output_dir) { save_file(output_dir + "species_count.csv"); }

// For species mode, save file directly in current directory.
void ReaxSpecies::save_file_to_current_dir() {
    save_file(file_path.substr(0, file_path.find_last_of(".")) + "_species_count.csv");
}
