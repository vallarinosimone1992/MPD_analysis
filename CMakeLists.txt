cmake_minimum_required(VERSION 3.16)
project(MPD_analysis LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

add_executable(reco src/reco_main.cpp)
target_include_directories(reco PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(reco PRIVATE ${ROOT_LIBRARIES})

add_executable(ana src/phys_main.cpp)
target_include_directories(ana PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(ana PRIVATE ${ROOT_LIBRARIES})
