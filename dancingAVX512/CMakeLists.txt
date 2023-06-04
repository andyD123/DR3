add_executable(dancingAVX512 dancingAVX512 AVX512Dance.cpp)

target_link_libraries(inverseCumNormalExample PUBLIC Vectorisation)

target_include_directories((dancingAVX512  PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )

        