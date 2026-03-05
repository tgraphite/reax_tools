cmake_minimum_required(VERSION 3.14)
project(reax_tools)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

file(GLOB_RECURSE SRC_FILES cpp/*.cpp)
list(APPEND SRC_FILES
    "${CMAKE_SOURCE_DIR}/cpp/fmt/format.cc"
    "${CMAKE_SOURCE_DIR}/cpp/fmt/os.cc"
)

add_executable(${PROJECT_NAME} ${SRC_FILES})


# Include directories for RDKIT and Boost
target_include_directories(${PROJECT_NAME} PRIVATE
    cpp
    cpp/fmt
    /opt/cpp/rdkit/include/rdkit
    /opt/cpp/boost_1_88_0/include/boost
)

# Link directories for RDKIT and Boost
target_link_directories(${PROJECT_NAME} PRIVATE
    /opt/cpp/rdkit/lib
    /opt/cpp/boost_1_88_0/lib
)

# Link RDKIT libraries (dynamic)
target_link_libraries(${PROJECT_NAME} PRIVATE
    RDKitGraphMol
    RDKitRDGeneral
    RDKitSmilesParse
    RDKitFileParsers
    RDKitRDGeometryLib
    RDKitRDStreams
    RDKitDistGeomHelpers
    RDKitMolDraw2D
)

