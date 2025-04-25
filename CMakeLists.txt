cmake_minimum_required(VERSION 3.21)

project(Mesh)

file(WRITE "build/config_file.txt" 
"# Config file for FEM solver

# Integration order
integration_order = 1

# Mesh file (located in meshes/)
mesh_file = square2d_4elt.msh

# Problem parameters
source = constant
boundary = zero

# Solution (optional)
sol = NONE
dx_sol = NONE
dy_sol = NONE

# Plot label
label = constant"
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Default build type" FORCE)
endif()

include(FetchContent)

# Fetch GoogleTest
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Fetch Eigen library (3.4.0)
set(EIGEN_BUILD_TESTS OFF)
FetchContent_Declare(
    Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0  
)
FetchContent_MakeAvailable(Eigen)

include_directories(include)

add_library(MeshUtilities "source/meshclass.cpp" "source/quadrature.cpp" "source/fem.cpp")
target_link_libraries(MeshUtilities PUBLIC Eigen3::Eigen)

add_executable(${PROJECT_NAME} "source/main.cpp")
target_link_libraries(${PROJECT_NAME} PRIVATE MeshUtilities Eigen3::Eigen)

#enable_testing()

#add_executable(hello_test "tests/hello_test.cpp")
#target_link_libraries(
#  hello_test
#  GTest::gtest_main
#)

#include(GoogleTest)
#gtest_discover_tests(hello_test)