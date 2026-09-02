add_executable(SmallBatchInference SmallBatchInference.cpp)

target_link_libraries(SmallBatchInference PUBLIC Vectorisation)

target_include_directories(SmallBatchInference PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}")

dr3_apply_isa(SmallBatchInference)

if(MSVC)
    target_compile_options(SmallBatchInference PRIVATE /EHsc /W4)
else()
    target_compile_options(SmallBatchInference PRIVATE -Wall -Wextra)
endif()

if(DR3_BUILD_TESTS)
    add_test(NAME SmallBatchInferenceSelfTest COMMAND SmallBatchInference --self-test)
    set_tests_properties(SmallBatchInferenceSelfTest PROPERTIES TIMEOUT 60)
endif()
