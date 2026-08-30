# SIMD portfolio pricing and Greeks

Standalone demonstration that prices a deterministic synthetic option
book with a scalar Black-Scholes reference and DR3 AVX2 (`VecD4D`).
It checks put-call parity, finite-difference Greeks, non-multiple-of-width
lengths, and rejects invalid vol/strike inputs.

This is an educational example, not a production pricer.

## Build and test

From the repository root, with tests enabled:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DDR3_BUILD_TESTS=ON -DDR3_BUILD_EXAMPLES=ON
cmake --build build --config Release --target portfolioGreeksExample
ctest --test-dir build -R portfolioGreeksSelfTest --output-on-failure
```

Or:

```
./build/portfolioGreeksExample/portfolioGreeksExample --self-test
```

## Benchmark

```
./build/portfolioGreeksExample/portfolioGreeksExample --size 4096
```

Printed times, throughput, and error maxima are measurements from the
local machine. They depend on processor, compiler, instruction set, and
build flags. The README does not claim a speedup.
