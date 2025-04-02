#include "getopt.h"

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// Modern-style getopt function.
// args: argc, argv, std::map<std::string, int> opt_nvals
// opt_nvals: e.g. {{"-f", 1}, {"--build", 0}, {"-s", 3}}
// return: std::map<std::string, std::vector<std::string>> opts_vals
// e.g. {{"-f", {"filename.txt"}}, {"--build", {}}, {"-s", {"1", "2", "3"}}}
std::map<std::string, std::vector<std::string>> neo_getopt(
    int argc, char *argv[], std::map<std::string, int> opts_nvals,
    const std::string &notif) {
    size_t cur_pos = 1;
    size_t cur_nvals = 0;
    std::string cur_arg;
    std::string cur_opt;
    std::map<std::string, std::vector<std::string>> opts_vals;

    while (cur_pos < argc) {
        cur_arg = argv[cur_pos];
        if (opts_nvals.find(cur_arg) != opts_nvals.end()) {
            cur_opt = cur_arg;
            cur_nvals = opts_nvals[cur_opt];
            opts_vals[cur_opt] = std::vector<std::string>();
        } else if (cur_nvals > 0) {
            opts_vals[cur_opt].push_back(cur_arg);
            cur_nvals--;
        } else {
            if (!notif.empty()) {
                std::cout << notif << std::endl;
            }
            throw std::runtime_error("Invalid options / parameters.");
        }
        cur_pos++;
    }

    return opts_vals;
}