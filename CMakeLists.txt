cmake_minimum_required(VERSION 3.17)

project(DR3)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


add_subdirectory(Vectorisation)
add_subdirectory(accumulateExample)
add_subdirectory(cumNormalExample)
add_subdirectory(inverseCumNormalExample)
