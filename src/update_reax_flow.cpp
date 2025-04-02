#include "reax_flow.h"
#include "universe.h"

void Universe::update_reax_flow(const std::shared_ptr<System> &prev_sys,
                                const std::shared_ptr<System> &curr_sys,
                                const int &curr_frame) {

    for (auto &prev_mol : prev_sys->molecules) {
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

        for (auto &curr_mol : curr_sys->molecules) {
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

            std::set_intersection(
                prev_mol->atom_ids.begin(), prev_mol->atom_ids.end(),
                curr_mol->atom_ids.begin(), curr_mol->atom_ids.end(),
                back_inserter(intersection));

            // Quick check: if intersection is empty, similarity is 0.
            if (intersection.empty())
                continue;

            // Quick check: if intersection size equals prev_mol size and
            // curr_mol size, then they are the same molecule.
            if (intersection.size() == prev_mol->atom_ids.size() &&
                intersection.size() == curr_mol->atom_ids.size())
                continue; // Same molecule, not reaction.

            std::set_union(prev_mol->atom_ids.begin(), prev_mol->atom_ids.end(),
                           curr_mol->atom_ids.begin(), curr_mol->atom_ids.end(),
                           back_inserter(union_set));

            float similarity =
                float(intersection.size()) / float(union_set.size());

            // If similarity exceeds the threshold, record and stop searching.
            if (similarity >= 0.5) {
                best_match = curr_mol;
                best_similarity = similarity;
                found_match = true;
                break; // If a good match is found, stop searching.
            }

            // If similarity exceeds the threshold, record and stop searching.
            else if (similarity > best_similarity) {
                best_match = curr_mol;
                best_similarity = similarity;
            }
        }

        // If found a match and not the same molecule (similarity between 0.5
        // and 1.0), record the reaction.
        if (best_match && best_similarity >= 0.5 && best_similarity < 1.0) {
            reax_flow->add_reaction(curr_frame, prev_mol, best_match);
        }
    }
}