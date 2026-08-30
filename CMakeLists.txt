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

function(dr3_apply_sse42 target)
    # VCL selects SSE4.2 Vec2d paths when INSTRSET >= 6. GCC/Clang get
    # that from -msse4.2. MSVC has no /arch:SSE4.2; /arch:AVX is the
    # lowest official arch that includes SSE4.2, plus INSTRSET=6.
    if(MSVC)
        target_compile_options(${target} PRIVATE /arch:AVX)
        target_compile_definitions(${target} PRIVATE INSTRSET=6)
    else()
        target_compile_options(${target} PRIVATE -msse4.2)
    endif()
endfunction()

function(dr3_apply_isa target)
    if(MSVC)
        if(DR3_ISA STREQUAL "AVX512")
            target_compile_options(${target} PUBLIC /arch:AVX512)
        elseif(DR3_ISA STREQUAL "SSE2")
            target_compile_options(${target} PUBLIC /arch:SSE2)
        elseif(DR3_ISA STREQUAL "SSE4.2" OR DR3_ISA STREQUAL "SSE42")
            target_compile_options(${target} PUBLIC /arch:AVX)
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

if(DR3_ENABLE_SANITIZERS)
    if(MSVC)
        message(WARNING "DR3_ENABLE_SANITIZERS is supported on GCC/Clang only")
    else()
        add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
        add_link_options(-fsanitize=address,undefined)
    endif()
endif()

if(DR3_BUILD_TESTS)
    include(CTest)
    enable_testing()
endif()

add_subdirectory(Vectorisation)

if(DR3_BUILD_EXAMPLES)
    add_subdirectory(accumulateExample)
    add_subdirectory(cumNormalExample)
    add_subdirectory(inverseCumNormalExample)
    add_subdirectory(dancingAVX512)
    add_subdirectory(lattice)
    add_subdirectory(GettingStarted)
    add_subdirectory(scratch)
    add_subdirectory(VariadicReduction)
    add_subdirectory(portfolioGreeksExample)
    add_subdirectory(portfolioStressExample)
endif()

if(DR3_BUILD_TESTS)
    add_subdirectory(VectorTest)
endif()

