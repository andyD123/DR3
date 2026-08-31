cmake_minimum_required(VERSION 3.17)

project(DR3 VERSION 0.2.0 LANGUAGES CXX)

# CMake package version files accept numeric versions only. Keep the release
# label separately so consumers and release tooling can distinguish this
# pre-release from the historical v01.3-alpha tag.
set(DR3_RELEASE_VERSION "0.2.0-alpha.1")

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(_DR3_BUILD_AUXILIARIES_DEFAULT OFF)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(_DR3_BUILD_AUXILIARIES_DEFAULT ON)
endif()

option(DR3_BUILD_TESTS "Build VectorTest and register it with CTest"
    ${_DR3_BUILD_AUXILIARIES_DEFAULT})
option(DR3_BUILD_EXAMPLES "Build example programs"
    ${_DR3_BUILD_AUXILIARIES_DEFAULT})
option(DR3_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer (GCC/Clang)" OFF)

unset(_DR3_BUILD_AUXILIARIES_DEFAULT)

# VectorTest uses the AVX2 DRC::VecD4D namespace. AVX512 may be compiled
# but is not the default runtime test ISA.
set(DR3_ISA "AVX2" CACHE STRING "Instruction set for CMake builds: SSE2, AVX2, or AVX512")
set_property(CACHE DR3_ISA PROPERTY STRINGS SSE2 AVX2 AVX512)

string(STRIP "${DR3_ISA}" DR3_ISA)
string(TOUPPER "${DR3_ISA}" DR3_ISA)
set(_DR3_SUPPORTED_ISAS SSE2 AVX2 AVX512)
if(NOT DR3_ISA IN_LIST _DR3_SUPPORTED_ISAS)
    string(JOIN ", " _DR3_SUPPORTED_ISAS_TEXT ${_DR3_SUPPORTED_ISAS})
    message(FATAL_ERROR
        "Unsupported DR3_ISA='${DR3_ISA}'. "
        "Supported values are: ${_DR3_SUPPORTED_ISAS_TEXT}")
endif()
set(DR3_ISA "${DR3_ISA}" CACHE STRING
    "Instruction set for CMake builds: SSE2, AVX2, or AVX512" FORCE)
set_property(CACHE DR3_ISA PROPERTY STRINGS SSE2 AVX2 AVX512)

function(dr3_apply_isa target)
    if(MSVC)
        if(DR3_ISA STREQUAL "AVX512")
            target_compile_options(${target} PUBLIC /arch:AVX512)
        elseif(DR3_ISA STREQUAL "SSE2")
            target_compile_options(${target} PUBLIC /arch:SSE2)
        elseif(DR3_ISA STREQUAL "AVX2")
            target_compile_options(${target} PUBLIC /arch:AVX2)
        else()
            message(FATAL_ERROR "Internal error: unhandled DR3_ISA='${DR3_ISA}'")
        endif()
    else()
        if(DR3_ISA STREQUAL "AVX512")
            target_compile_options(${target} PUBLIC -mavx512f -mfma)
        elseif(DR3_ISA STREQUAL "SSE2")
            target_compile_options(${target} PUBLIC -msse2)
        elseif(DR3_ISA STREQUAL "AVX2")
            target_compile_options(${target} PUBLIC -mavx2 -mfma)
        else()
            message(FATAL_ERROR "Internal error: unhandled DR3_ISA='${DR3_ISA}'")
        endif()
    endif()
endfunction()

function(dr3_enable_sanitizers target)
    if(DR3_ENABLE_SANITIZERS AND NOT MSVC)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined
            -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined)
    endif()
endfunction()

if(DR3_ENABLE_SANITIZERS AND MSVC)
    message(WARNING "DR3_ENABLE_SANITIZERS is supported on GCC/Clang only")
endif()

add_subdirectory(Vectorisation)
dr3_enable_sanitizers(Vectorisation)

if(DR3_BUILD_EXAMPLES)
    add_subdirectory(accumulateExample)
    add_subdirectory(cumNormalExample)
    add_subdirectory(inverseCumNormalExample)
    add_subdirectory(dancingAVX512)
    add_subdirectory(lattice)
    add_subdirectory(GettingStarted)
    add_subdirectory(scratch)
    add_subdirectory(VariadicReduction)

    foreach(dr3_example_target IN ITEMS
            accumulateExample
            cumNormalExample
            inverseCumNormalExample
            dancingAVX512
            lattice
            GettingStarted
            ScratchExample
            VariadicExample)
        dr3_enable_sanitizers(${dr3_example_target})
    endforeach()
endif()

if(DR3_BUILD_TESTS)
    include(CTest)
    enable_testing()
    add_subdirectory(VectorTest)
    dr3_enable_sanitizers(VectorTest)
endif()
set(DR3_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/DR3")

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/DR3Config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/DR3Config.cmake"
    INSTALL_DESTINATION "${DR3_INSTALL_CMAKEDIR}")

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/DR3ConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMinorVersion)

install(EXPORT DR3Targets
    FILE DR3Targets.cmake
    NAMESPACE DR3::
    DESTINATION "${DR3_INSTALL_CMAKEDIR}")

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/DR3Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/DR3ConfigVersion.cmake"
    DESTINATION "${DR3_INSTALL_CMAKEDIR}")

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/DR3")
