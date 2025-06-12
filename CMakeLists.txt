set(PROJECT_NAME reax_tools)
cmake_minimum_required(VERSION 3.10)

project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Force static linking
set(BUILD_SHARED_LIBS OFF)
add_subdirectory(external/fmt)

file(GLOB_RECURSE SRC_FILES "src/*.cpp" "src/*.c")
file(GLOB_RECURSE HEADER_FILES "include/*.h" "external/fmt/include/*.h")

add_executable(${PROJECT_NAME} ${SRC_FILES})

target_include_directories(${PROJECT_NAME} PRIVATE 
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(${PROJECT_NAME} PRIVATE fmt::fmt)

# Add compiler flags
target_compile_options(${PROJECT_NAME} PRIVATE -O3 -march=native -mtune=native -mavx2 -mfma -g -std=c++17 -static)

# Add linker flags for static linking of standard libraries
target_link_options(${PROJECT_NAME} PRIVATE -static)
