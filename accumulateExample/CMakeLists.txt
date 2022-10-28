add_executable(accumulateExample accumulate_example.cpp)

target_link_libraries(accumulateExample PUBLIC Vectorisation)

target_include_directories(accumulateExample PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )