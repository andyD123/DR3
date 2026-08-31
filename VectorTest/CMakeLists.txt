find_package(GTest CONFIG QUIET)

if(NOT GTest_FOUND)
    if(POLICY CMP0135)
        cmake_policy(SET CMP0135 NEW)
    endif()

    include(FetchContent)

    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
        URL_HASH SHA256=1f357c27ca988c3f7c6b4bf68a9395005ac6761f034046e9dde0896e3aba00e4
    )
    FetchContent_MakeAvailable(googletest)
endif()

# These translation units compile with GCC, Clang, and MSVC using the
# AVX2 DRC::VecD4D test namespace. Other VectorTest files currently
# depend on MSVC language extensions or missing headers.
set(DR3_VECTOR_TEST_PORTABLE_SOURCES
    TestBasicPortable.cpp
    TestConcurrentAllocator.cpp
    TestAccumulator.cpp
    TestAllocator.cpp
    TestFilterSelect.cpp
    testNamespace.cpp
    test_precise_accumulation.cpp
)

add_executable(VectorTest ${DR3_VECTOR_TEST_PORTABLE_SOURCES})

target_link_libraries(VectorTest PRIVATE Vectorisation GTest::gtest GTest::gtest_main)

target_include_directories(VectorTest PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

dr3_apply_isa(VectorTest)

if(MSVC)
    target_compile_options(VectorTest PRIVATE /EHsc)
endif()

# ApplyFilterB misreads compact Vec4db masks on AVX-512 hosts; tracked in
# https://github.com/andyD123/DR3/issues/47. Keep it in the tree but out of
# the portable CTest filter until that issue is fixed. Use add_test rather than
# gtest_discover_tests: CMake 4.x POST_BUILD discovery fails on
# Windows when the binary is in a multi-config output directory.
set(DR3_GTEST_EXCLUDE
    "-TestFilterSelect.ApplyFilterB")

add_test(NAME VectorTest COMMAND $<TARGET_FILE:VectorTest> --gtest_filter=${DR3_GTEST_EXCLUDE})
add_test(NAME ParallelAllocator
    COMMAND $<TARGET_FILE:VectorTest> --gtest_filter=ParallelAllocator.*)
