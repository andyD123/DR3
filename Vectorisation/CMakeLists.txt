add_library(Vectorisation STATIC
            VecX/alloc_policy.cpp
            VecX/vec.cpp)
add_library(DR3::Vectorisation ALIAS Vectorisation)

dr3_apply_isa(Vectorisation)

set_target_properties(Vectorisation PROPERTIES
    POSITION_INDEPENDENT_CODE ON)

target_compile_features(Vectorisation PUBLIC cxx_std_17)

target_include_directories(Vectorisation PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

install(TARGETS Vectorisation
    EXPORT DR3Targets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")

install(DIRECTORY
    "${CMAKE_CURRENT_SOURCE_DIR}/VecX"
    "${CMAKE_CURRENT_SOURCE_DIR}/VCL"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp")

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/VCL/LICENSE"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/DR3"
    RENAME VCL-LICENSE)
