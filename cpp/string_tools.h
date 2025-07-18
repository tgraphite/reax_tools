#pragma once
#include <iostream>
#include <map>
#include <string>

#include "defines.h"

// Convert formula to element-number map.
// e.g. C1H2 -> {{"C": 1}, {"H": 2}}
std::map<std::string, int> parse_formula(const std::string& formula);

// Rename and normalize formula by element order,
// e.g. order = C,H,O,N
// C1O2N1H3 -> CH3NO2
std::string rename_formula(const std::string& formula, const std::vector<std::string>& order = ELEMENT_DISPLAY_ORDER);

std::vector<std::string> split_by_space(const std::string& str);
std::vector<std::string> split(const std::string& str, const std::string& delim);

bool can_convert_to_int(const std::string& str);
bool starts_with(const std::string& str, const std::string& prefix);
bool ends_with(const std::string& str, const std::string& suffix);

// This template function accepts a container of pointer (e.g.
// vector<shared_ptr<Actual-type>>) prints the ptr->info() for the first N
// elements. if two_way = true, prints first and last N prt->info().
template <typename T> void print_info(const T& container, const size_t& num, const bool& two_way = true) {
    size_t size = container.size();

    if (size == 0)
        return;

    if (two_way) {
        if (size < num * 2) {
            for (size_t i = 0; i < size; i++) {
                std::cout << (container[i])->info();
            }
            return;
        }
        for (size_t ileft = 0; ileft < num; ileft++) {
            std::cout << (container[ileft])->info();
        }
        for (size_t iright = size - num; iright < size; iright++) {
            std::cout << (container[iright])->info();
        }
    } else {
        if (size < num) {
            for (size_t i = 0; i < size; i++) {
                std::cout << (container[i])->info();
                return;
            }
        }
        for (size_t i = 0; i < num; i++) {
            std::cout << (container[i])->info();
        }
    }
    std::cout << "..." << std::endl;
}

// operator<< for std::vector<std::string>
std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& vs);