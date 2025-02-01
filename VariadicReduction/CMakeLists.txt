add_executable(VariadicExample VariadicReduction.cpp)

target_link_libraries(VariadicExample PUBLIC Vectorisation)

target_include_directories(VariadicExample PUBLIC
        "${PROJECT_BINARY_DIR}"
)
