cmake_minimum_required(VERSION 3.10)

set(CMAKE_BUILD_TYPE Debug)

project(jsondotrulo VERSION 0.01)

configure_file(src/c_make_header.h.in src/c_make_header.h)

add_subdirectory(third_party)

add_executable(jsondotrulo
                src/main.cc
                src/graph.cc)

target_include_directories(jsondotrulo PUBLIC
                           "${PROJECT_BINARY_DIR}/src"
                           )

target_link_libraries(jsondotrulo PRIVATE nlohmann_json::nlohmann_json)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED true)
