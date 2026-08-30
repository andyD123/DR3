add_executable(portfolioStressExample portfolio_stress.cpp)

target_link_libraries(portfolioStressExample PUBLIC Vectorisation)

target_include_directories(portfolioStressExample PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_SOURCE_DIR}/portfolioGreeksExample")

if(DR3_BUILD_TESTS)
    add_test(NAME portfolioStressSelfTest COMMAND portfolioStressExample --self-test)
    set_tests_properties(portfolioStressSelfTest PROPERTIES TIMEOUT 60)
endif()
