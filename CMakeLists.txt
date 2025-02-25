cmake_minimum_required(VERSION 3.21)

project(Mesh)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Default build type" FORCE)
endif()

#find_package(Eigen3 REQUIRED)
#if (Eigen3_FOUND)
#    message(STATUS "Eigen found")
#endif()

include_directories(include)

add_library(MeshUtilities "source/meshclass.cpp")
target_link_libraries(MeshUtilities PUBLIC) # Eigen3::Eigen)

# main executable
add_executable(${PROJECT_NAME} "source/main.cpp")
target_link_libraries(${PROJECT_NAME} PRIVATE MeshUtilities) # Eigen3::Eigen)