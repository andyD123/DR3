cmake_minimum_required(VERSION 3.17)

project(DR3SubdirectoryConsumer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

if(NOT DEFINED DR3_ROOT)
    message(FATAL_ERROR "Configure with -DDR3_ROOT=/path/to/DR3")
endif()

add_subdirectory("${DR3_ROOT}" dr3 EXCLUDE_FROM_ALL)

if(DR3_BUILD_TESTS OR DR3_BUILD_EXAMPLES)
    message(FATAL_ERROR
        "DR3 tests and examples must default to OFF under add_subdirectory")
endif()

foreach(_unexpected_target IN ITEMS VectorTest accumulateExample lattice)
    if(TARGET "${_unexpected_target}")
        message(FATAL_ERROR
            "Unexpected auxiliary DR3 target: ${_unexpected_target}")
    endif()
endforeach()

add_executable(dr3_subdirectory_consumer main.cpp)
target_link_libraries(dr3_subdirectory_consumer
    PRIVATE DR3::Vectorisation)

enable_testing()
add_test(NAME DR3SubdirectoryConsumer
    COMMAND $<TARGET_FILE:dr3_subdirectory_consumer>)
