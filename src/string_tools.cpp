#include "string_tools.h"

#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// GPT-4 generated. Chemical formula to map<element: number>
std::map<std::string, int> parse_formula(const std::string &formula) {
    std::map<std::string, int> elements_nums;
    std::string element;
    std::string number;

    for (size_t i = 0; i < formula.size(); ++i) {
        if (isalpha(formula[i])) {
            if (!element.empty() && !number.empty()) {
                elements_nums[element] += std::stoi(number);
                number.clear();
            } else if (!element.empty()) {
                elements_nums[element] += 1;  // when no number in formula (e.g.
                                              // CO), make default = 1.
            }

            element = formula[i];  // new element

            // for symbols with two chars (e.g. Fe)
            if (islower(formula[i + 1])) {
                element += formula[++i];
            }
        } else if (isdigit(formula[i])) {
            number += formula[i];
        }
    }

    if (!element.empty()) {
        if (number.empty()) number = "1";
        elements_nums[element] += std::stoi(number);  // the last element.
    }

    return elements_nums;
}

// format formula by order.
// default order: see sorted_elements in constants.h
std::string rename_formula(const std::string &formula,
                           const std::vector<std::string> &order) {
    std::map<std::string, int> elements_nums = parse_formula(formula);
    std::string result;

    for (const auto &elem : order) {
        auto it = elements_nums.find(elem);
        if (it != elements_nums.end()) {
            result += it->first;
            if (it->second != 1)
                result += std::to_string(it->second);
            else
                result += "1";
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

std::vector<std::string> split(const std::string &str,
                               const std::string &delim) {
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