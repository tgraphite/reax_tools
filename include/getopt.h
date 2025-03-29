#pragma once

#include <map>
#include <vector>
#include <string>

std::map<std::string, std::vector<std::string>> neo_getopt(int argc, char* argv[], std::map<std::string, int> opts_nvals, const std::string& notif = "");