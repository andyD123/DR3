add_executable(cumNormalExample cumNormalExample.cpp)

target_link_libraries(cumNormalExample PUBLIC Vectorisation)

target_include_directories(cumNormalExample PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )