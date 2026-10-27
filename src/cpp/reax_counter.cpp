#include "reax_counter.h"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <set>

#include "fmt/format.h"
#include "string_tools.h"

void SpeciesCounter::import_frame_species(int frame_id, const std::vector<SpeciesRecord>& frame_species) {
    if (frame_id < 1) return;
    nframes = std::max(nframes, frame_id);

    for (const auto& species : frame_species) {
        auto it = records.find(species.id);
        if (it == records.end()) {
            records.emplace(species.id, species);
        } else if (!species.smiles.empty() && !it->second.smiles.empty() &&
                   species.smiles != it->second.smiles) {
            fmt::print(stderr,
                       "Warning: structure ID collision for id={} ({} vs {})\n",
                       species.id, it->second.smiles, species.smiles);
        }

        auto& values = counts[species.id];
        if (values.size() < static_cast<std::size_t>(frame_id)) {
            values.resize(frame_id, 0);
        }
        values[frame_id - 1] += 1;
    }
}

void SpeciesCounter::analyze_frame_species() {
    for (auto& [id, values] : counts) {
        if (values.size() < static_cast<std::size_t>(nframes)) {
            values.resize(nframes, 0);
        }
    }
}

void SpeciesCounter::save_species_count() {
    FILE* file = create_file("species_count.csv");
    if (file == nullptr) return;

    fmt::print(file, "frame");
    for (const auto& [id, record] : records) {
        fmt::print(file, ",{}", id);
    }
    fmt::print(file, "\n");

    for (int frame = 0; frame < nframes; ++frame) {
        fmt::print(file, "{}", frame + 1);
        for (const auto& [id, record] : records) {
            auto it = counts.find(id);
            int value = 0;
            if (it != counts.end() && static_cast<std::size_t>(frame) < it->second.size()) {
                value = it->second[frame];
            }
            fmt::print(file, ",{}", value);
        }
        fmt::print(file, "\n");
    }
    fclose(file);
}

void SpeciesCounter::save_molecules() {
    FILE* file = create_file("molecules.csv");
    if (file == nullptr) return;

    fmt::print(file, "id,formula,smiles\n");
    for (const auto& [id, record] : records) {
        fmt::print(file, "{},{},{}\n", id, csv_escape(record.formula), csv_escape(record.smiles));
    }
    fclose(file);
}

std::size_t SpeciesCounter::formula_count() const {
    std::set<std::string> formulas;
    for (const auto& [id, record] : records) formulas.insert(record.formula);
    return formulas.size();
}

void SpeciesCounter::brief_report() {
    if (records.empty()) {
        fmt::print("=== Species Report ===\n(no species)\n");
        return;
    }

    std::vector<std::pair<double, const SpeciesRecord*>> ranking;
    ranking.reserve(records.size());
    for (const auto& [id, record] : records) {
        double total = 0.0;
        auto it = counts.find(id);
        if (it != counts.end()) {
            total = std::accumulate(it->second.begin(), it->second.end(), 0.0);
        }
        ranking.emplace_back(total, &record);
    }
    std::sort(ranking.begin(), ranking.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    const std::size_t limit = std::min<std::size_t>(ranking.size(), 20);
    fmt::print("=== Species Report ===\n");
    fmt::print("{:<20s}{:>10s}  {}\n", "formula", "average", "smiles");
    for (std::size_t i = 0; i < limit; ++i) {
        const auto* record = ranking[i].second;
        fmt::print("{:<20s}{:>10.2f}  {}\n", record->formula,
                   ranking[i].first / std::max(nframes, 1), record->smiles);
    }
}
