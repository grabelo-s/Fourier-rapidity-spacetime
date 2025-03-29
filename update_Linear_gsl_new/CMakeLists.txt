cmake_minimum_required(VERSION 3.10)

project(Main VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/includes/>
)

file(GLOB sources ./src/*.cpp)
file(GLOB headers ./includes/*.hpp)

configure_file(./input/xv_Q2_CT14_positive_valence.dat ./xv_Q2_CT14_positive_valence.dat COPYONLY)
configure_file(./input/xv_Q2_CT14_positive_valence_ud_sea_ud.dat ./xv_Q2_CT14_positive_valence_ud_sea_ud.dat COPYONLY)
configure_file(./input/gaussian_F.dat ./gaussian_F.dat COPYONLY)

find_package(GSL REQUIRED)
add_executable(main main.cpp ${sources} ${headers})

target_link_libraries(main GSL::gsl)
