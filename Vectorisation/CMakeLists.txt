add_library(Vectorisation STATIC
            VecX/alloc_policy.cpp
            VecX/vec.cpp)

if(MSVC)
    dr3_apply_isa(Vectorisation)
else()
    # Keep historical GCC/Clang behaviour so AVX512 examples still compile
    # on machines where -march=native includes those instructions.
    if(DR3_ISA STREQUAL "AVX512")
        target_compile_options(Vectorisation PUBLIC -mavx512f -mfma)
    elseif(DR3_ISA STREQUAL "SSE2")
        target_compile_options(Vectorisation PUBLIC -msse2)
    else()
        target_compile_options(Vectorisation PUBLIC -march=native)
    endif()
endif()

target_include_directories(Vectorisation PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}")

            