set(PROJECT_NAME reax_tools)
cmake_minimum_required(VERSION 3.10)

project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)

set(BUILD_SHARED_LIBS ON)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")

# Set Boost path
set(BOOST_ROOT "/opt/cpp/boost_1_88_0")
set(BOOST_INCLUDEDIR "/opt/cpp/boost_1_88_0/include")
set(BOOST_LIBRARYDIR "/opt/cpp/boost_1_88_0/lib")

# Set RDKit path
set(RDKit_DIR "/opt/cpp/rdkit/lib/cmake/rdkit")
set(RDKit_INCLUDE_DIR "/opt/cpp/rdkit/include/rdkit")
set(RDKit_LIBRARY_DIR "/opt/cpp/rdkit/lib")

find_package(RDKit REQUIRED)
include_directories(${RDKit_INCLUDE_DIR})

file(GLOB_RECURSE SRC_FILES "src/*.cpp" "src/*.c")
file(GLOB_RECURSE HEADER_FILES "include/*.h")

add_subdirectory(external/fmt)

add_executable(${PROJECT_NAME} ${SRC_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME}
    PRIVATE
    fmt::fmt
    ${RDKit_LIBRARIES}
    -L${RDKit_LIBRARY_DIR}
    RDKitRDGeneral
    RDKitGraphMol
    RDKitSmilesParse
    RDKitFileParsers
    RDKitMolDraw2D
    RDKitDistGeometry
    RDKitDistGeomHelpers
)

# Add compiler flags
target_compile_options(${PROJECT_NAME} PRIVATE -O3 -march=native -mtune=native -mavx2 -mfma -g -std=c++17)

# Add linker flags for static linking of standard libraries
target_link_options(${PROJECT_NAME} PRIVATE -static-libstdc++ -static-libgcc)
