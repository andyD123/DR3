cmake_minimum_required(VERSION 3.17)

project(DR3)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(DR3_BUILD_TESTS "Build VectorTest and register it with CTest" ON)
option(DR3_BUILD_EXAMPLES "Build example programs" ON)
option(DR3_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer (GCC/Clang)" OFF)

# VectorTest uses the AVX2 DRC::VecD4D namespace. AVX512 may be compiled
# but is not the default runtime test ISA.
set(DR3_ISA "AVX2" CACHE STRING "Instruction set for CMake builds: SSE2, SSE4.2, AVX2, or AVX512")
set_property(CACHE DR3_ISA PROPERTY STRINGS SSE2 SSE4.2 AVX2 AVX512)

function(dr3_apply_isa target)
    if(MSVC)
        if(DR3_ISA STREQUAL "AVX512")
            target_compile_options(${target} PUBLIC /arch:AVX512)
        elseif(DR3_ISA STREQUAL "SSE2")
            target_compile_options(${target} PUBLIC /arch:SSE2)
        elseif(DR3_ISA STREQUAL "SSE4.2" OR DR3_ISA STREQUAL "SSE42")
            # Keep MSVC's normal x64 SSE2 code-generation baseline. VCL uses
            # explicit SSE4.2 intrinsics when INSTRSET is 6, so /arch:AVX is
            # neither necessary nor compatible with an SSE4.2 minimum.
            target_compile_definitions(${target} PUBLIC INSTRSET=6)
        else()
            target_compile_options(${target} PUBLIC /arch:AVX2)
        endif()
    else()
        if(DR3_ISA STREQUAL "AVX512")
            target_compile_options(${target} PUBLIC -mavx512f -mfma)
        elseif(DR3_ISA STREQUAL "SSE2")
            target_compile_options(${target} PUBLIC -msse2)
        elseif(DR3_ISA STREQUAL "SSE4.2" OR DR3_ISA STREQUAL "SSE42")
            target_compile_options(${target} PUBLIC -msse4.2)
        else()
            target_compile_options(${target} PUBLIC -mavx2 -mfma)
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
