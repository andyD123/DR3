# Building DR3

## Requirements and instruction sets

DR3 requires CMake 3.17, a C++17 compiler, and an x86 or x64 target. The
required build and package paths are exercised with MSVC, GCC, and Clang on
Windows and Ubuntu.

`DR3_ISA` controls the compiler baseline for DR3 and every target that links
its public headers. Accepted values are `SSE2`, `AVX2`, and `AVX512`; other
values stop configuration with an error. The default is `AVX2`, which maps to
`-mavx2 -mfma` with GCC and Clang and `/arch:AVX2` with MSVC. DR3 does not use
`-march=native`, so builds do not silently depend on the build machine.

This option selects compile flags, not runtime dispatch or a public vector
type. Source code continues to select an existing namespace such as
`DRC::VecD4D`. AVX512 is opt-in, and an AVX512 binary must not be executed
unless the host CPU supports it.

## Build from source

```sh
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DDR3_BUILD_TESTS=ON \
  -DDR3_BUILD_EXAMPLES=ON \
  -DDR3_ISA=AVX2
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Available options are:

- `DR3_BUILD_EXAMPLES=ON|OFF` (default `ON` standalone, `OFF` as a subdirectory)
- `DR3_BUILD_TESTS=ON|OFF` (default `ON` standalone, `OFF` as a subdirectory)
- `DR3_ISA=SSE2|AVX2|AVX512` (default `AVX2`)
- `DR3_ENABLE_SANITIZERS=ON|OFF` (default `OFF`; GCC/Clang ASan and UBSan)

Enabling tests fetches the pinned GoogleTest dependency. A core-only build has
no test dependency:

```sh
cmake -S . -B build-core \
  -DCMAKE_BUILD_TYPE=Release \
  -DDR3_BUILD_TESTS=OFF \
  -DDR3_BUILD_EXAMPLES=OFF
cmake --build build-core --config Release
```

## Install and use the CMake package

Build and install one instruction-set variant per prefix. The installed target
exports the same compile baseline to consumers of DR3's inline headers. The
static core is position-independent so it can also be linked into shared
consumer libraries.

```sh
cmake -S . -B build-package \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PWD/dr3-install" \
  -DDR3_BUILD_TESTS=OFF \
  -DDR3_BUILD_EXAMPLES=OFF \
  -DDR3_ISA=AVX2
cmake --build build-package --config Release
cmake --install build-package --config Release
```

A consuming project can then use the versioned package and namespaced target:

```cmake
find_package(DR3 0.2 CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE DR3::Vectorisation)
```

```cpp
#include <VecX/dr3.h>
```

`<VecX/dr3.h>` is the supported installed umbrella header. The package also
installs the headers needed transitively by that entry point, but standalone
inclusion of every internal `VecX` header is not a compatibility promise.

DR3 is still pre-1.0. This package carries the release label
`0.2.0-alpha.1`; CMake's numeric package version is `0.2.0`. Its version file
uses same-minor compatibility, so `0.2.x` packages can satisfy a
`find_package(DR3 0.2 ...)` request; compatibility across different pre-1.0
minor versions is not implied.

The package exposes `DR3_BUILT_ISA` as informational metadata. Consumers must
not change it; install another DR3 build to a separate prefix when a different
baseline is required.

The checked-in smoke test links the installed archive into a shared consumer
library and runs a small executable against it. This verifies position-
independent code, headers, exported compile requirements, package version, and
the `find_package` entry point:

```sh
cmake -S tests/package_consumer -B build-package-consumer \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PWD/dr3-install" \
  -DDR3_EXPECTED_ISA=AVX2
cmake --build build-package-consumer --config Release
ctest --test-dir build-package-consumer -C Release --output-on-failure
```

For an in-tree dependency, the same stable target name is available without an
installation:

```cmake
add_subdirectory(path/to/DR3)
target_link_libraries(my_target PRIVATE DR3::Vectorisation)
```

Tests and examples default to `OFF` in this form so DR3 does not fetch
GoogleTest or add auxiliary executables to the parent build. Set either option
before `add_subdirectory` when those targets are wanted explicitly. The
checked-in `tests/subdirectory_consumer` fixture verifies these defaults.

## Test scope

The portable CTest suite uses the AVX2 `DRC::VecD4D` namespace and therefore
requires AVX2 and FMA at runtime. The pairwise transform-reduce tail regressions
are included in the portable suite. `TestFilterSelect.ApplyFilterB` remains
excluded because it fails outside the original MSVC Release configuration.
