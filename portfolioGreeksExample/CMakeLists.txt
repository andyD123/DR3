add_executable(portfolioGreeksExample portfolio_greeks.cpp)

target_link_libraries(portfolioGreeksExample PUBLIC Vectorisation)

target_include_directories(portfolioGreeksExample PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}")

if(DR3_BUILD_TESTS)
    add_test(NAME portfolioGreeksSelfTest COMMAND portfolioGreeksExample --self-test)
    set_tests_properties(portfolioGreeksSelfTest PROPERTIES TIMEOUT 60)
endif()
