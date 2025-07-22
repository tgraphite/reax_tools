#include "argparser.h"
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    // Initialize parser
    ArgParser parser = init_argparser();

    // Parse command line arguments
    try {
        parser.parse_args(argc, argv);
        parser.operate_for_all();
    } catch (const std::exception& e) {
        std::cerr << "ArgParser error: " << e.what() << std::endl;
        return 1;
    }

    // Print all global variables for verification
    std::cout << "===== ArgParser Global Variables =====" << std::endl;
    std::cout << "INPUT_FILE: " << INPUT_FILE << std::endl;
    std::cout << "INPUT_ELEMENT_TYPES: ";
    for (const auto& elem : INPUT_ELEMENT_TYPES) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
    std::cout << "RVDW_FACTOR: " << RVDW_FACTOR << std::endl;
    std::cout << "NUM_THREADS: " << NUM_THREADS << std::endl;
    std::cout << "MAX_REACTIONS: " << MAX_REACTIONS << std::endl;
    std::cout << "FLAG_DUMP_STRUCTURE: " << FLAG_DUMP_STRUCTURE << std::endl;
    std::cout << "DUMP_STEPS: " << DUMP_STEPS << std::endl;
    std::cout << "FLAG_MARK_RING_ATOMS: " << FLAG_MARK_RING_ATOMS << std::endl;
    std::cout << "FLAG_NO_RINGS: " << FLAG_NO_RINGS << std::endl;
    std::cout << "FLAG_NO_REACTIONS: " << FLAG_NO_REACTIONS << std::endl;
    std::cout << "FLAG_NO_REDUCE_REACTIONS: " << FLAG_NO_REDUCE_REACTIONS << std::endl;
    std::cout << "MERGE_TARGET: " << MERGE_TARGET << std::endl;
    std::cout << "MERGE_RANGES: ";

    for (const auto& v : MERGE_RANGES) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    std::cout << "FLAG_RESCALE_MERGE_COUNT: " << FLAG_RESCALE_MERGE_COUNT << std::endl;
    std::cout << "ELEMENT_DISPLAY_ORDER: " << std::endl;

    for (const auto& elem : ELEMENT_DISPLAY_ORDER) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;

    std::cout << "ELEMENT_ATOMIC_RADII: " << std::endl;
    for (const auto& [element, val] : ELEMENT_ATOMIC_RADII) {
        std::cout << element << "=" << val << " ";
    }
    std::cout << std::endl;

    std::cout << "ELEMENT_VALENCES: " << std::endl;
    for (const auto& [element, val] : ELEMENT_MAX_VALENCIES) {
        std::cout << element << "=" << val << " ";
    }
    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    return 0;
}
