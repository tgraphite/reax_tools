#pragma once

#include <string_tools.h>

#include <filesystem>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

class ArgParser {
  private:
    struct ArgOption {
        std::string name;                                  // 选项名
        std::string short_name;                            // 短选项名
        std::string description;                           // 描述
        std::string default_value;                         // 默认值
        bool required;                                     // 是否必需
        bool is_flag;                                      // 是否为标志（无需参数值）
        std::string group;                                 // 分组
        std::string value_type;                            // 值类型描述
        std::function<bool(const std::string&)> validator; // 值验证函数
    };

    std::string program_name;
    std::string program_description;
    std::map<std::string, ArgOption> options;
    std::map<std::string, std::string> parsed_values;
    std::vector<std::string> positional_args;
    std::map<std::string, std::vector<std::string>> groups;
    std::string example_usage;

  public:
    ArgParser(const std::string& name, const std::string& description)
        : program_name(name), program_description(description) {}

    void add_example(const std::string& example) { example_usage = example; }

    void add_argument(const std::string& name, const std::string& short_name, const std::string& description,
                      const std::string& group = "General", const std::string& default_value = "",
                      bool required = false, bool is_flag = false, const std::string& value_type = "",
                      std::function<bool(const std::string&)> validator = nullptr) {
        ArgOption option;
        option.name = name;
        option.short_name = short_name;
        option.description = description;
        option.default_value = default_value;
        option.required = required;
        option.is_flag = is_flag;
        option.group = group;
        option.value_type = value_type;
        option.validator = validator;

        options[name] = option;
        if (!short_name.empty()) {
            options[short_name] = option;
        }

        if (!group.empty()) {
            groups[group].push_back(name);
        }
    }

    bool parse_args(int argc, char* argv[]) {
        if (argc <= 1) {
            print_help();
            return false;
        }

        std::vector<std::string> args(argv + 1, argv + argc);
        size_t i = 0;

        while (i < args.size()) {
            std::string arg = args[i];

            if (arg == "-h" || arg == "--help") {
                print_help();
                return false;
            }

            if (arg.substr(0, 2) == "--" || arg.substr(0, 1) == "-") {
                std::string option_name = arg;

                if (options.find(option_name) == options.end()) {
                    std::cerr << "Error: Invalid option " << option_name << std::endl;
                    return false;
                }

                const ArgOption& option = options[option_name];

                if (option.is_flag) {
                    parsed_values[option.name] = "true";
                    i++;
                } else {
                    if (i + 1 >= args.size()) {
                        std::cerr << "Error: option " << option_name << " need values." << std::endl;
                        return false;
                    }

                    std::string value = args[i + 1];

                    // 检验参数值
                    if (option.validator && !option.validator(value)) {
                        std::cerr << "Error: option " << option_name << " got invalid value '" << value << std::endl;
                        return false;
                    }

                    parsed_values[option.name] = value;
                    i += 2;
                }
            } else {
                positional_args.push_back(arg);
                i++;
            }
        }

        // 检查必需参数
        for (const auto& pair : options) {
            const ArgOption& option = pair.second;
            if (option.name != pair.first)
                continue;
            if (option.required &&
                (parsed_values.find(option.name) == parsed_values.end() || parsed_values[option.name].empty())) {
                std::cerr << "Error: Mandatory arguments missing! " << option.name << std::endl;
                return false;
            }
        }

        return true;
    }

    // 打印帮助信息
    void print_help() const {
        const int DESC_WIDTH = 120;
        const int OPTION_WIDTH = 30;

        std::cout << program_description << std::endl << std::endl;
        std::cout << "Usages: " << program_name << " [options]" << std::endl << std::endl;

        // 按分组打印选项
        for (const auto& group_pair : groups) {
            std::cout << group_pair.first << ":" << std::endl;

            for (const auto& option_name : group_pair.second) {
                const ArgOption& option = options.at(option_name);

                if (option.name != option_name)
                    continue;

                std::stringstream option_str;
                if (!option.short_name.empty()) {
                    option_str << option.short_name << ", ";
                }
                option_str << option.name;

                if (!option.is_flag) {
                    if (!option.value_type.empty()) {
                        option_str << " <" << option.value_type << "> ";
                    } else {
                        option_str << " VALUE";
                    }
                }

                std::string option_display = option_str.str();

                std::cout << "  " << std::left << std::setw(OPTION_WIDTH) << option_display;

                // 处理说明文字的换行
                std::string desc = option.description;
                if (!option.default_value.empty()) {
                    desc += " (default: " + option.default_value + ")";
                }

                if (desc.length() <= DESC_WIDTH - OPTION_WIDTH - 2) {
                    std::cout << desc << std::endl;
                } else {
                    size_t pos = 0;
                    bool first_line = true;

                    while (pos < desc.length()) {
                        size_t len = std::min(DESC_WIDTH - (first_line ? OPTION_WIDTH + 2 : 8),
                                              static_cast<int>(desc.length() - pos));
                        size_t space_pos = desc.rfind(' ', pos + len);

                        if (space_pos > pos && space_pos < pos + len) {
                            len = space_pos - pos;
                        }

                        if (!first_line) {
                            std::cout << "        ";
                        }

                        std::cout << desc.substr(pos, len) << std::endl;

                        pos += len;
                        if (pos < desc.length() && desc[pos] == ' ')
                            pos++;
                        first_line = false;
                    }
                }
            }

            std::cout << std::endl;
        }

        if (!example_usage.empty()) {
            std::cout << "Examples:" << std::endl;
            std::cout << "  " << example_usage << std::endl << std::endl;
        }
    }

    template <typename T> T get(const std::string& name) const {
        if (parsed_values.find(name) == parsed_values.end()) {
            const ArgOption& option = options.at(name);
            if (option.is_flag) {
                return get_default_for_flag<T>();
            }
            if (!option.default_value.empty()) {
                return convert<T>(option.default_value);
            }
            return T();
        }

        return convert<T>(parsed_values.at(name));
    }

    // 为标志类型选项提供适当的默认值
    template <typename T> T get_default_for_flag() const {
        return T(); // 默认构造
    }

    bool has_flag(const std::string& name) const {
        return parsed_values.find(name) != parsed_values.end() && parsed_values.at(name) == "true";
    }
    bool has_option(const std::string& name) const { return parsed_values.find(name) != parsed_values.end(); }
    const std::vector<std::string>& get_positional_args() const { return positional_args; }

  private:
    template <typename T> T convert(const std::string& value) const {
        std::istringstream ss(value);
        T result;
        ss >> result;
        return result;
    }
};

// 外部定义模板特化
template <> inline bool ArgParser::get_default_for_flag<bool>() const { return false; }

template <> inline std::string ArgParser::get_default_for_flag<std::string>() const { return "false"; }

template <> inline std::string ArgParser::convert<std::string>(const std::string& value) const { return value; }

template <> inline bool ArgParser::convert<bool>(const std::string& value) const {
    return value == "true" || value == "1" || value == "yes";
}

template <>
inline std::vector<std::string> ArgParser::convert<std::vector<std::string>>(const std::string& value) const {
    return split(value, ",");
}

template <> inline std::vector<int> ArgParser::convert<std::vector<int>>(const std::string& value) const {
    std::vector<std::string> parts = split(value, ",");
    std::vector<int> result;
    for (const auto& part : parts) {
        result.push_back(std::stoi(part));
    }
    return result;
}

bool validate_float_positive(const std::string& value) {
    try {
        float val = std::stof(value);
        return val > 0;
    } catch (...) {
        return false;
    }
}

bool validate_int_positive(const std::string& value) {
    try {
        int val = std::stoi(value);
        return val > 0;
    } catch (...) {
        return false;
    }
}