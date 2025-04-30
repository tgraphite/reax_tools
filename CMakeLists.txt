set(PROJECT_NAME reax_tools)
cmake_minimum_required(VERSION 3.8)

project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
set(BUILD_SHARED_LIBS OFF)

file(GLOB_RECURSE SRC_FILES "src/*.cpp" "src/*.c")
file(GLOB_RECURSE HEADER_FILES "include/*.h")

add_subdirectory(external/fmt)

add_executable(${PROJECT_NAME} ${SRC_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} PRIVATE fmt::fmt) 

if (UNIX AND NOT APPLE)
    target_compile_options(${PROJECT_NAME} PRIVATE -O3 -march=native -mtune=native -mavx2 -mfma -g -static-libstdc++ -static-libgcc)
    target_compile_definitions(${PROJECT_NAME} PRIVATE _LINUX)
    target_link_options(${PROJECT_NAME} PRIVATE -static)
elseif (WIN32)
    target_compile_definitions(${PROJECT_NAME} PRIVATE _WIN32 )
 #   target_compile_options(${PROJECT_NAME} PRIVATE /MT)
endif()
