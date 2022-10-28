add_executable(inverseCumNormalExample inverseCumNormalExample.cpp cdfNormalInverse.cpp)

target_link_libraries(inverseCumNormalExample PUBLIC Vectorisation)

target_include_directories(inverseCumNormalExample PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )

        