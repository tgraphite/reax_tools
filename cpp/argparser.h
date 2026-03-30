/// @file argparser.h
/// @brief Provides ArgParser, a very strict command-line argument parser for
/// scientific computing applications.

#pragma once

#include <fmt/core.h>

#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

extern int MAX_RING_SIZE;
extern int MIN_RING_SIZE;
extern int MAX_NEIGH;

extern int NUM_THREADS;
extern int MAX_REACTIONS;
extern int MAX_KEY_MOLECULES;
extern int DUMP_STEPS;

extern std::string INPUT_FILE;
extern std::string OUTPUT_DIR;
extern std::vector<std::string> INPUT_ELEMENT_TYPES;
extern float RVDW_FACTOR;

extern bool FLAG_DUMP_STRUCTURE;
extern bool FLAG_MARK_RING_ATOMS;

extern bool FLAG_NO_REDUCE_REACTIONS;
extern bool FLAG_NO_RINGS;
extern bool FLAG_NO_REACTIONS;

// Reaction tracking flags (ChemTrayzer-like) - ON by default
extern bool FLAG_TRACK_REACTIONS;
extern bool FLAG_NO_TRACK_REACTIONS;
extern int STABLE_TIME_FRAMES;
extern float TIMESTEP_FS;
extern int SAMPLING_FREQ;

extern bool FLAG_RESCALE_MERGE_COUNT;
extern std::string MERGE_TARGET;
extern std::vector<int> MERGE_RANGES;

extern std::vector<std::string> ALL_ELEMENTS;
extern std::map<std::string, int> ELEMENT_TO_INDEX;
extern std::map<std::string, float> ELEMENT_ATOMIC_RADII;
extern std::map<std::string, float> ELEMENT_MASS;
extern std::map<std::string, int> ELEMENT_MAX_VALENCIES;
extern std::vector<std::string> ELEMENT_DISPLAY_ORDER;
extern std::vector<unsigned int> PRIME_NUMBERS;
extern std::vector<unsigned int> BIGGER_PRIME_NUMBERS;

/**
 * @class ArgParserException
 * @brief Exception type for all errors thrown by ArgParser.
 *
 * This exception is thrown whenever an error occurs during argument parsing or
 * operation execution.
 */
class ArgParserException : public std::runtime_error {
   public:
    /**
     * @brief Construct a new ArgParserException with a message.
     * @param msg The error message.
     */
    explicit ArgParserException(const std::string& msg) : std::runtime_error(msg) {}
};

/**
 * @class ArgParser
 * @brief Strict command-line argument parser for scientific computing
 * applications.
 *
 * ArgParser enforces that every option must have an associated operation and
 * that all arguments are strictly validated. It is designed for use cases where
 * input correctness is critical, such as scientific computing.
 *
 * Usage example:
 * @code
 * ArgParser parser("myprog", "Program description");
 * parser.add_argument("--input", "-i", "Input file", "File", false,
 * [](std::string val){ ... }); parser.parse_args(argc, argv);
 * parser.operate_for_all();
 * @endcode
 */
class ArgParser {
   private:
    /**
     * @struct Option
     * @brief Represents a single command-line option and its associated metadata
     * and operation(s).
     *
     * Each Option must have at least one operation defined. Operations can be for
     * single value, multiple values, or flag.
     */
    struct Option {
        std::string short_name;                  // Short option name (e.g., -i)
        std::string long_name;                   // Long option name (e.g., --input)
        std::string desc;                        // Description of the option
        std::string desc_group;                  // Group name for help display
        std::vector<std::string> parsed_values;  // Values parsed from command line
        bool is_flag;                            // True if this option is a flag (no value expected)
        bool used = false;
        bool closed = false;

        std::function<bool(std::vector<std::string>)> operation;  // Operation for single value

        // Default constructor
        Option() : short_name(""), long_name(""), desc(""), desc_group(""), is_flag(false), operation(nullptr) {}

        /**
         * @brief Construct a new Option object.
         * @param short_name_ Short option name
         * @param long_name_ Long option name
         * @param desc_ Description
         * @param desc_group_ Group name
         * @param is_flag_ True if flag
         * @param operation Operation for multiple values (std::vector<std::string>)
         * @throws std::logic_error if no operation is provided
         */
        Option(std::string short_name_, std::string long_name_, std::string desc_, std::string desc_group_,
               bool is_flag_, std::function<bool(std::vector<std::string>)> operation_ = nullptr)
            : short_name(short_name_),
              long_name(long_name_),
              desc(desc_),
              desc_group(desc_group_),
              is_flag(is_flag_),
              operation(operation_) {}

        /**
         * @brief Execute the associated operation for this option.
         * @return true if operation succeeded, false otherwise
         * @throws ArgParserException if no parameter is set or no valid operation
         * is found
         */
        bool operate() {
            if (is_flag) return operation(parsed_values);

            if (parsed_values.size() == 0) {
                throw ArgParserException(
                    fmt::format("Error when processing options: Parameter not set for option {}", long_name));
            } else if (operation == nullptr) {
                throw ArgParserException(
                    fmt::format("Error when processing options: Can not find operation for {}", long_name));
            } else {
                return operation(parsed_values);
            }
        }
    };

    std::string program_name;                                     // Name of the program
    std::string program_desc;                                     // Description of the program
    std::map<std::string, Option> options;                        // Map of long option name to Option
    std::map<std::string, std::vector<std::string>> desc_groups;  // Map of group name to option names

   public:
    /**
     * @brief Construct a new ArgParser object.
     * @param _program_name Name of the program
     * @param _program_desc Description of the program
     */
    ArgParser(std::string _program_name, std::string _program_desc)
        : program_name(_program_name), program_desc(_program_desc) {}

    // -------------------- Developer phase --------------------

    /**
     * @brief Add a new argument (option) to the parser.
     *
     * At least one operation must be provided, otherwise a logic_error will be
     * thrown.
     *
     * @param long_name Long option name (e.g., --input)
     * @param short_name Short option name (e.g., -i)
     * @param desc Description of the option
     * @param desc_group Group name for help display
     * @param is_flag True if this option is a flag (no value expected)
     * @param operation Operation for multiple values (std::vector<std::string>)
     * @throws std::logic_error if no operation is provided
     */
    void add_argument(std::string long_name, std::string short_name, std::string desc, std::string desc_group,
                      bool is_flag, std::function<bool(std::vector<std::string>)> operation) {
        Option option(short_name, long_name, desc, desc_group, is_flag, operation);
        options[long_name] = option;
        desc_groups[desc_group].push_back(long_name);
    }

    /**
     * @brief Get the canonical long name for an option given any valid name.
     * @param arg The option name (long or short)
     * @return The canonical long name
     * @throws ArgParserException if the argument is invalid
     */
    std::string get_option_long_name(std::string arg) {
        for (const auto& [long_name, option] : options) {
            if (option.long_name == arg || option.short_name == arg) {
                return option.long_name;
            }
        }

        throw ArgParserException(fmt::format("Error when processing arguments: invalid argument {}", arg));
    }

    /**
     * @brief Print help information for the program and its options.
     *
     * @note Will invoke exit(0) to force user retry after reading help
     * information.
     */
    void print_help() {
        fmt::print("Usage: {} [options]\n", program_name);
        fmt::print("{}\n\n", program_desc);

        for (const auto& [group, option_names] : desc_groups) {
            fmt::print("[{}]\n", group);
            for (const auto& long_name : option_names) {
                const auto& opt = options.at(long_name);
                // Format:   -s, --long   Description
                fmt::print("    {:<6}, {:<16}: {}\n", opt.short_name, opt.long_name, opt.desc);
            }
            fmt::print("\n");
        }
        exit(0);
    }

// -------------------- User phase --------------------

/**
 * @brief Parse command-line arguments.
 *
 * This function parses the command-line arguments and fills the options
 * accordingly.
 *
 * @param argc Argument count
 * @param argv Argument vector
 * @return true if parsing succeeded, false otherwise
 * @throws ArgParserException on any parsing error
 */
#ifndef WASM_MODE
    bool parse_args(int argc, char* argv[])
#else
    bool parse_args(int argc, const char* argv[])
#endif
    {
        if (argc <= 1) {
            print_help();
            return false;
        }

        // Convert argv to vector<string>
        std::vector<std::string> args(argv + 1, argv + argc);

        size_t tmp_index;
        std::string curr_option = "";

        for (size_t index = 0; index < args.size(); ++index) {
            std::string curr_arg = args[index];

            // An option name
            if (curr_arg.substr(0, 2) == "--" || curr_arg.substr(0, 1) == "-") {
                if (!curr_arg.empty() && options.find(curr_arg) != options.end() && !options[curr_option].is_flag) {
                    options[curr_option].closed = true;
                }
                curr_option = get_option_long_name(curr_arg);
                if (options[curr_option].is_flag) {
                    options[curr_option].used = true;
                    options[curr_option].closed = true;
                }
            } else {  // A parameter
                if (curr_option.empty()) {
                    throw ArgParserException(
                        fmt::format("Error when processing arguments: Can not find "
                                    "option for parameter {}",
                                    curr_option));
                }
                if (options[curr_option].closed) {
                    throw ArgParserException(
                        fmt::format("Error when processing arguments: Set option {} multiple times", curr_option));
                }
                if (!options[curr_option].is_flag) {
                    options[curr_option].parsed_values.push_back(curr_arg);
                    options[curr_option].used = true;
                } else {
                    throw ArgParserException(
                        fmt::format("Error when processing arguments: Can not set "
                                    "parameter for flag {}",
                                    curr_option));
                }
            }
        }

        return true;
    }

    /**
     * @brief Execute the operation for all options.
     *
     * This function will call the associated operation for each option. If any
     * operation fails, an exception is thrown.
     *
     * @throws ArgParserException if any operation fails
     */
    void operate_for_all() {
        for (auto& [long_name, option] : options) {
            if (!option.used) {
                continue;
            }

            bool has_successed = option.operate();

            if (!has_successed) {
                throw ArgParserException(
                    fmt::format("Error when processing option {}, check format and "
                                "numbers of parameters",
                                option.long_name));
            }
        }
    }
};

ArgParser init_argparser();
