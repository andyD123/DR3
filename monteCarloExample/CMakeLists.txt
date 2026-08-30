add_executable(monteCarloExample monte_carlo.cpp)

target_link_libraries(monteCarloExample PUBLIC Vectorisation)

target_include_directories(monteCarloExample PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_SOURCE_DIR}/portfolioGreeksExample")

if(DR3_BUILD_TESTS)
    add_test(NAME monteCarloSelfTest COMMAND monteCarloExample --self-test)
    set_tests_properties(monteCarloSelfTest PROPERTIES TIMEOUT 120)
endif()
