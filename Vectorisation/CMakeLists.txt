add_library(Vectorisation STATIC
            VecX/alloc_policy.cpp
            VecX/vec.cpp)

if (MSVC)    
    # add_compile_options(/W4 /WX)
else()
    target_compile_options(Vectorisation PUBLIC "-march=native")

    #  or -mavx/-mavx2/-mavx512f (and -march= options that imply them with relevant tuning settings)
    #target_compile_options(Vectorisation PUBLIC "--std=c++17")
    #target_compile_options(Vectorisation PUBLIC "-mavx2")
    #target_compile_options(Vectorisation PUBLIC "-mfma")
endif()

            