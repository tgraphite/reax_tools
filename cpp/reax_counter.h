#pragma once
#include <fmt/core.h>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "argparser.h"
#include "string_tools.h"

template <typename T> class Counter {
  private:
    std::map<std::string, std::vector<T>> entry_frame_values;
    int nframes = 0;

  public:
    Counter() = default;

    // void import_frame_values(const int& frame_id, const std::map<std::string, T>& entry_values);
    // void save_file(const std::string& output_path);

    void import_frame_values(const int& frame_id, const std::map<std::string, T>& entry_values) {
        size_t index = frame_id - 1; // frame_id starts from 1

        for (const auto& [entry, value] : entry_values) {
            if (entry_frame_values.find(entry) == entry_frame_values.end()) {
                std::vector<T> new_frame_values;
                new_frame_values.resize(index + 1, T());
                entry_frame_values[entry] = new_frame_values;
            }

            if (entry_frame_values[entry].size() < index + 1) {
                entry_frame_values[entry].resize(index + 1, T());
            }

            entry_frame_values[entry][index] = value;
        }
    }
    void save_file(const std::string& basename) {
        if (entry_frame_values.empty())
            return;

        FILE* fp = create_file(basename);

        // Write header
        bool first = true;
        for (const auto& pair : entry_frame_values) {
            if (!first)
                fmt::print(fp, ",");
            fmt::print(fp, "{}", pair.first);
            first = false;
        }
        fmt::print(fp, "\n");

        // Find max frame count
        size_t max_frames = 0;
        for (const auto& pair : entry_frame_values) {
            if (pair.second.size() > max_frames)
                max_frames = pair.second.size();
        }

        // Write each frame
        for (size_t i = 0; i < max_frames; ++i) {
            first = true;
            for (const auto& pair : entry_frame_values) {
                if (!first)
                    fmt::print(fp, ",");
                if (i < pair.second.size())
                    fmt::print(fp, "{}", pair.second[i]);
                else
                    fmt::print(fp, "{}", T());
                first = false;
            }
            fmt::print(fp, "\n");
        }

        fclose(fp);
    }
};

class SpeciesCounter {
  private:
    std::string file_path;
    std::string save_path;
    std::ifstream file;
    std::ofstream save;
    std::vector<std::string> formulas;

    std::map<std::string, std::vector<float>> formulas_nums;
    std::vector<std::vector<std::string>> all_frame_formulas;
    int nframes = 0;

  public:
    // When read from Lammps reax/c/bonds file.
    SpeciesCounter(const std::string& file_path);

    // For import data from class Universe.
    SpeciesCounter();
    ~SpeciesCounter();

    void get_formulas();
    void get_nums();
    void import_frame_formulas(int& frame_id, const std::vector<std::string>& frame_formulas);
    void analyze_frame_formulas();

    // The --order option
    void rename_all_formulas(const std::vector<std::string>& order);

    // The -me option
    void merge_formulas(const std::vector<std::string>& formulas, const std::string& new_formula);
    void merge_by_element(const std::string& target_element, const std::vector<int>& ranges, bool recalc = false);

    // The -rc option
    void scale_formula(const std::string& formula, const float& k);
    void rescale_all_by_element(const std::string& target_element);

    void show_nums();
    void brief_report();

    void save_file();
};
