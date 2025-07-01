cmake_minimum_required(VERSION 3.10)
project(reax_tools)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXE_LINKER_FLAGS "-static")  # 强制全静态

find_package(RDKit REQUIRED)
find_package(Boost REQUIRED COMPONENTS serialization iostreams system filesystem thread regex date_time chrono atomic timer program_options random math_c99 math_tr1)

# 自动查找所有 .cpp 文件（不包括 .cc）
file(GLOB_RECURSE PROJECT_CPP_FILES
    "${CMAKE_SOURCE_DIR}/cpp/*.cpp"
)
# 排除 .cc 文件
list(FILTER PROJECT_CPP_FILES EXCLUDE REGEX "\\.cc$")

# 显式指定 .cc 文件
set(PROJECT_CC_FILES
    cpp/fmt/format.cc
    cpp/fmt/os.cc
)

# 合并所有源文件
set(ALL_SRC_FILES
    ${PROJECT_CPP_FILES}
    ${PROJECT_CC_FILES}
    cpp/universe_parallel.cc
    cpp/local_main.cc
    cpp/rdkit_utils.cc
)

# 打印调试信息
message(STATUS "All source files: ${ALL_SRC_FILES}")

add_executable(reax_tools ${ALL_SRC_FILES})

target_include_directories(reax_tools PRIVATE cpp cpp/fmt)
target_compile_definitions(reax_tools PRIVATE ENABLE_RDKIT)

target_link_libraries(reax_tools
    PRIVATE
        RDKit::GraphMol
        RDKit::RDGeneral
        RDKit::SmilesParse
        RDKit::FileParsers
        RDKit::RDGeometryLib
        RDKit::RDStreams
        RDKit::DistGeomHelpers
        RDKit::MolDraw2D
        Boost::serialization
        Boost::iostreams
        Boost::system
        Boost::filesystem
        Boost::thread
        Boost::regex
        Boost::date_time
        Boost::chrono
        Boost::atomic
        Boost::timer
        Boost::program_options
        Boost::random
        Boost::math_c99
        Boost::math_tr1
        freetype
)
