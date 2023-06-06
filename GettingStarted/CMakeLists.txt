add_executable(GettingStarted GettingStarted.cpp)

target_link_libraries(GettingStarted PUBLIC Vectorisation)

target_include_directories(GettingStarted PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )