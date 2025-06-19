#pragma once

#include <map>
#include <string>
#include <vector>

std::map<std::string, std::vector<std::string>> neo_getopt(
    int argc, char* argv[], std::map<std::string, int> opts_nvals,
    const std::string& notif = "");