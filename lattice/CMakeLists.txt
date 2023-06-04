add_executable(lattice lattice.cpp)

target_link_libraries(lattice PUBLIC Vectorisation)

target_include_directories(lattice PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )

        