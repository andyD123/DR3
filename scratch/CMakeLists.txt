add_executable(ScratchExample scratch.cpp)

target_link_libraries(ScratchExample PUBLIC Vectorisation)

target_include_directories(ScratchExample PUBLIC
        "${PROJECT_BINARY_DIR}"
)