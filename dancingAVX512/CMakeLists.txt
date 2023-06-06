add_executable(dancingAVX512 dancingAVX512.cpp AVX512Dance.cpp)

target_link_libraries(dancingAVX512 PUBLIC Vectorisation)

target_include_directories(dancingAVX512  PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )

        