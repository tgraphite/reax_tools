#include "string_tools.h"

#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "fmt/format.h"

// GPT-4 generated. Chemical formula to map<element: number>
std::map<std::string, int> parse_formula(const std::string &formula) {
    std::map<std::string, int> elements_nums;
    std::string element;
    std::string number;
    size_t i = 0;

    while (i < formula.size()) {
        // Start of a new element
        if (isalpha(formula[i])) {
            // Process previous element if exists
            if (!element.empty()) {
                int count = number.empty() ? 1 : std::stoi(number);
                elements_nums[element] += count;
                number.clear();
            }

            // Start new element
            element = formula[i];

            // Check for two-character element (e.g., Fe, Cl)
            if (i + 1 < formula.size() && islower(formula[i + 1])) {
                element += formula[i + 1];
                i++;
            }
        }
        // Accumulate number
        else if (isdigit(formula[i])) {
            number += formula[i];
        }

        i++;
    }

    // Process the last element
    if (!element.empty()) {
        int count = number.empty() ? 1 : std::stoi(number);
        elements_nums[element] += count;
    }

    return elements_nums;
}

// format formula by order.
// default order: see sorted_elements in constants.h
std::string rename_formula(const std::string &formula, const std::vector<std::string> &order) {
    std::map<std::string, int> elements_nums = parse_formula(formula);
    std::string result;

    if (starts_with(formula, "grp_")) {
        return formula;
    }

    // First add elements in the specified order
    for (const auto &elem : order) {
        auto it = elements_nums.find(elem);
        if (it != elements_nums.end()) {
            result += it->first;
            if (it->second > 1) {
                result += std::to_string(it->second);
            }
            elements_nums.erase(it);  // Remove processed element
        }
    }

    // Then add any remaining elements in alphabetical order
    for (const auto &pair : elements_nums) {
        result += pair.first;
        if (pair.second > 1) {
            result += std::to_string(pair.second);
        }
    }

    return result;
}

// Other function.
std::vector<std::string> split_by_space(const std::string &str) {
    std::istringstream iss(str);
    std::vector<std::string> tokens;
    std::string token;

    tokens.reserve(20);

    while (iss >> token) {
        tokens.push_back(token);
    }

    return tokens;
}

std::vector<std::string> split(const std::string &str, const std::string &delim) {
    std::vector<std::string> tokens;
    size_t prev = 0, pos = 0;
    while (pos < str.length() && prev < str.length()) {
        pos = str.find(delim, prev);
        if (pos == std::string::npos) pos = str.length();
        std::string token = str.substr(prev, pos - prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    return tokens;
}

/// Using try-catch exception will make too much debug pollution, and cost a
/// little bit more time.
// This function is used to check if the first column of xyz file is a number or
// element symbol.
bool can_convert_to_int(const std::string &str) {
    for (size_t i = 0; i < str.size(); i++) {
        if (!isdigit(str[i])) return false;
    }
    return true;
}

std::ostream &operator<<(std::ostream &os, const std::vector<std::string> &vs) {
    for (const auto &v : vs) {
        os << v << ", ";
    }
    return os;
}

// old implementation of starts_with
bool starts_with(const std::string &str, const std::string &prefix) {
    return str.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), str.begin());
}

// old implementation of ends_with
bool ends_with(const std::string &str, const std::string &suffix) {
    return str.size() >= suffix.size() && std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}
