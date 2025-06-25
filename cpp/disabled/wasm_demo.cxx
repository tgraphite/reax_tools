/**
 * This cpp file is for learning and testing WebAssembly embedding.
 * Mock a classic command-line program behavior: read input files, do something, save output file.
 * And the raw main() function is replaced by an interface compatible with javascript.
 * Author: Graphite, Date: 2025.06.20
 */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// extern "C" int cpp_main(int argc, const char **argv);

// /*
// int main()
// {
//     return 0;
// } // 空main
// */

// 实际业务逻辑（模拟你的命令行程序）
extern "C" int cpp_main(int argc, const char **argv)
{
    std::string input_filename;
    std::string s_param;
    std::string t_param;

    // 1. 简易的命令行参数解析
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "-f")
        {
            if (i + 1 < argc)
            {
                input_filename = argv[++i];
            }
            else
            {
                std::cerr << "错误: -f 参数需要一个文件名" << std::endl;
                return 1;
            }
        }
        else if (arg == "-s")
        {
            if (i + 1 < argc)
            {
                s_param = argv[++i];
            }
            else
            {
                std::cerr << "错误: -s 参数需要一个值" << std::endl;
                return 1;
            }
        }
        else if (arg == "-t")
        {
            if (i + 1 < argc)
            {
                t_param = argv[++i];
            }
            else
            {
                std::cerr << "错误: -t 参数需要一个值" << std::endl;
                return 1;
            }
        }
    }

    // 2. 检查必要参数
    if (input_filename.empty())
    {
        std::cerr << "用法: " << argv[0] << " -f <文件名> [-s <值>] [-t <值>]" << std::endl;
        return 1;
    }

    // 3. 处理 -s 和 -t 参数
    if (!s_param.empty())
    {
        std::cout << "提示: 收到参数 -s, 值为: " << s_param << std::endl;
    }
    if (!t_param.empty())
    {
        std::cout << "提示: 收到参数 -t, 值为: " << t_param << std::endl;
    }

    // 4. 读取输入文件
    std::cerr << "尝试打开输入文件: " << input_filename << std::endl;
    std::ifstream input(input_filename);
    if (!input)
    {
        std::cerr << "错误: 打开输入文件失败" << std::endl;
        return 1;
    }

    // 5. 处理数据（示例：转大写）
    std::string content((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    for (char &c : content)
        if (c >= 'a' && c <= 'z')
            c -= 32;

    // 6. 写入输出文件
    std::string output_filename = "/output/result.txt";
    std::ofstream output(output_filename);
    output << content;
    output.close();
    std::cout << "已处理 " << content.length() << " 字符，并写入到 " << output_filename << std::endl;

    // 7. 创建另一个日志文件，以测试多文件 zip 下载
    std::string log_filename = "/output/log.txt";
    std::ofstream log_file(log_filename);
    log_file << "处理摘要:\n";
    log_file << "输入文件: " << input_filename << "\n";
    log_file << "输出文件: " << output_filename << "\n";
    if (!s_param.empty())
    {
        log_file << "S-参数: " << s_param << "\n";
    }
    if (!t_param.empty())
    {
        log_file << "T-参数: " << t_param << "\n";
    }
    log_file.close();
    std::cout << "已创建日志文件 " << log_filename << std::endl;

    return 0;
}
