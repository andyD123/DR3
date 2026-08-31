# DR3

DR3 is a C++17 SIMD container and algorithm layer. Elementwise operations,
reductions, and predicates are written as generic lambdas; filters turn those
predicates into owning selections of matching values.

The current API chooses scalar type and register width through a namespace such
as `DRC::VecD4D`. `DR3_ISA` selects the compiler baseline for the complete
target. Both choices are compile-time decisions: DR3 does not currently provide
runtime CPU dispatch.

## Why this layer exists

The vendored VCL2 backend supplies SIMD register types. DR3 adds an owned,
padded `Vec`, owning filtered selections, and container-shaped transform,
filter, and reduction algorithms. It is intended for fixed-ISA numeric code
that benefits from expressing the operation once as a lambda. It is not a
tensor framework or a runtime-dispatch system.

## A small example

This example deliberately has five inputs, so the AVX2 four-double path also
exercises a non-vector-width tail.

```cpp
#include <VecX/dr3.h>

#include <vector>

int main()
{
    using Vec = DRC::VecD4D::VecXX;

    Vec values(std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0});
    auto square = [](auto x) { return x * x; };
    auto squared = transform(square, values);

    auto atLeastTen = [](auto x) { return x >= Vec::INS(10.0); };
    auto selected = filter(atLeastTen, squared);

    return selected.size() == 2 &&
                   selected[0] == 16.0 && selected[1] == 25.0
        ? 0
        : 1;
}
```

`selected` owns a snapshot of the matching values; later changes to `squared`
are not reflected in it.

## ISA model

The namespace identifies the scalar and register shape used by the source code:

| Namespace | Scalar | Lanes | Intended native ISA |
| --- | --- | ---: | --- |
| `DRC::VecD2D` | `double` | 2 | SSE2, 128-bit |
| `DRC::VecF4F` | `float` | 4 | SSE2, 128-bit |
| `DRC::VecD4D` | `double` | 4 | AVX2, 256-bit |
| `DRC::VecF8F` | `float` | 8 | AVX2, 256-bit |
| `DRC::VecD8D` | `double` | 8 | AVX512, 512-bit |
| `DRC::VecF16F` | `float` | 16 | AVX512, 512-bit |

`DR3_ISA=SSE2|AVX2|AVX512` controls compiler flags. The default is AVX2 plus
FMA (`-mavx2 -mfma` with GCC/Clang or `/arch:AVX2` with MSVC); builds never
fall back to `-march=native`. Wider VCL register shapes may be emulated at a
lower baseline, but there is no runtime selection between implementations.

Keep the namespace and compiler baseline consistent across translation units
that exchange DR3 types. AVX512 is opt-in: compile-only checks are safe on
other hosts, but never execute an AVX512 binary without confirming CPU support.

## Build from source

```sh
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DDR3_BUILD_TESTS=ON \
  -DDR3_BUILD_EXAMPLES=OFF \
  -DDR3_ISA=AVX2
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

## Install and consume

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

```cmake
find_package(DR3 0.2 CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE DR3::Vectorisation)
```

`<VecX/dr3.h>` is the supported installed umbrella header. This package is the
`0.2.0-alpha.1` pre-release; its numeric CMake version is `0.2.0` and it
provides same-minor (`0.2.x`) compatibility.

## Current contract and caveats

- **Tails:** core portable tests cover logical lengths that are not multiples
  of the register width, including pairwise transform-reduce tails. Only
  `[0, size())` is public data; padded elements are an implementation detail.
  Experimental/internal algorithms do not yet have the same coverage.
- **Alignment:** pool allocation sizes are rounded to cache-line-sized element
  counts, and each pooled block starts at a 64-byte boundary.
- **Filtered results:** the `VecView` result owns copied values and their source
  indices. It does not borrow the input, so source mutation is not reflected in
  an existing result.
- **Concurrency:** shared mutation is unsupported. The allocator registry on
  this branch is unsynchronized, so concurrent operations on independently
  owned `Vec` values are not yet guaranteed safe. Allocator synchronization is
  tracked separately in [PR #42](https://github.com/andyD123/DR3/pull/42).
- **Dispatch:** ISA selection is compile-time only. There is no CPUID-based
  runtime fallback.

See [Building DR3](docs/Build.md) for package details, supported options, and
the compiler/CI matrix. DR3 is licensed under the [Apache License 2.0](LICENSE).
