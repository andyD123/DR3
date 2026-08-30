# Supported platforms

OS: Windows, Ubuntu
Compiler: Visual Studio, gcc, clang
CPU: AVX2 is the default CMake instruction set (`DR3_ISA=AVX2`). The
VectorTest suite uses the `DRC::VecD4D` namespace and therefore needs
AVX2 at runtime. AVX512 can be compiled with `-DDR3_ISA=AVX512` but
those binaries should not be executed unless the CPU supports AVX512.

# Build commands

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Optional CMake flags:

- `-DDR3_BUILD_EXAMPLES=ON|OFF` (default ON)
- `-DDR3_BUILD_TESTS=ON|OFF` (default ON)
- `-DDR3_ISA=SSE2|AVX2|AVX512` (default AVX2)
- `-DDR3_ENABLE_SANITIZERS=ON` (GCC/Clang AddressSanitizer and
  UndefinedBehaviorSanitizer)

# Tests

The CTest target uses the AVX2 `DRC::VecD4D` namespace. The portable
suite is the VectorTest translation units that compile with GCC, Clang,
and MSVC, plus `TestBasicPortable.cpp`. Additional VectorTest files are
compiled on MSVC only; `TestCurve.cpp` is omitted because
`ExampleVectors/curve.h` is not in the repository.

Two pairwise transform-reduce cases currently crash under GCC/Clang and
are excluded from the portable CTest filter.

```
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DDR3_BUILD_TESTS=ON

cmake --build build --config Release

ctest --test-dir build -C Release --output-on-failure
```

Example-only configure path:

```
cmake -S . -B build -DDR3_BUILD_TESTS=OFF
```

