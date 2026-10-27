#pragma once
#include <fmt/core.h>
#include <cstdint>
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
        fmt::print(fp, "frame,");

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
            fmt::print(fp, "{},", i + 1); // Frame number starts from 1

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

struct SpeciesRecord {
    std::uint64_t id = 0;
    std::string formula;
    std::string smiles;  // empty for formula-level fallback
};

class SpeciesCounter {
private:
    std::map<std::uint64_t, SpeciesRecord> records;
    std::map<std::uint64_t, std::vector<int>> counts;
    int nframes = 0;

public:
    SpeciesCounter() = default;
    ~SpeciesCounter() = default;

    void import_frame_species(int frame_id, const std::vector<SpeciesRecord>& frame_species);
    void analyze_frame_species();
    void save_species_count();
    void save_molecules();
    void brief_report();

    std::size_t species_count() const { return records.size(); }
    std::size_t formula_count() const;
};
