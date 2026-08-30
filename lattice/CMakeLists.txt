add_library(dr3_lattice STATIC
    american_trinomial.cpp
    barrier_trinomial.cpp
    european_binomial.cpp
    european_trinomial.cpp
    tree_pricer_detail.cpp
    tridiagonal.cpp
)

target_link_libraries(dr3_lattice PUBLIC Vectorisation)
target_include_directories(dr3_lattice PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
dr3_apply_isa(dr3_lattice)

if(DR3_BUILD_EXAMPLES)
    add_executable(latticeExample lattice_example.cpp)
    target_link_libraries(latticeExample PRIVATE dr3_lattice)
endif()

if(DR3_BUILD_TESTS)
    add_executable(LatticeTest
        tests/TestPdeGrid.cpp
        tests/TestTreePricers.cpp
        tests/TestTridiagonal.cpp
    )
    target_link_libraries(LatticeTest PRIVATE dr3_lattice GTest::gtest GTest::gtest_main)
    dr3_apply_isa(LatticeTest)
    if(MSVC)
        target_compile_options(LatticeTest PRIVATE /EHsc)
    endif()
    add_test(NAME LatticeTest COMMAND $<TARGET_FILE:LatticeTest>)
    add_test(NAME Grid1D COMMAND $<TARGET_FILE:LatticeTest> --gtest_filter=Grid1D.*)
    add_test(NAME Tridiagonal COMMAND $<TARGET_FILE:LatticeTest> --gtest_filter=Tridiagonal.*)
endif()
