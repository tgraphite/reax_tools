cmake_minimum_required(VERSION 3.14)
project(reax_tools)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
endif()

option(REAX_ENABLE_THREADS "Enable multi-threaded frame processing" ON)
option(REAX_FORCE_STDCXXFS "Force linking libstdc++fs (for old GCC compatibility testing)" OFF)

file(GLOB_RECURSE SRC_FILES src/cpp/*.cpp)
list(APPEND SRC_FILES
    "${CMAKE_SOURCE_DIR}/src/cpp/fmt/format.cc"
    "${CMAKE_SOURCE_DIR}/src/cpp/fmt/os.cc"
)

add_executable(reax_tools_core ${SRC_FILES})

target_include_directories(reax_tools_core PRIVATE
    src/cpp
    src/cpp/fmt
)

if(REAX_ENABLE_THREADS)
    find_package(Threads REQUIRED)
    target_compile_definitions(reax_tools_core PRIVATE REAX_ENABLE_THREADS=1)
    target_link_libraries(reax_tools_core PRIVATE Threads::Threads)
endif()

# Test-only SMILES targets.  The production core does not depend on them.
add_executable(smiles_selftest test/cpp/smiles_selftest.cpp src/cpp/smiles.cpp)
target_include_directories(smiles_selftest PRIVATE src/cpp)

add_executable(smiles_cli test/cpp/smiles_cli.cpp src/cpp/smiles.cpp)
target_include_directories(smiles_cli PRIVATE src/cpp)

# GCC 8.x keeps std::filesystem in a separate library; GCC 9+ does not.
if(REAX_FORCE_STDCXXFS OR (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9))
    target_link_libraries(reax_tools_core PRIVATE stdc++fs)
endif()
