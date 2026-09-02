# Small-batch MLP inference

This example loads one committed float32 `5 -> 7 -> 3` model and runs a
two-sample batch through:

```
LayerNorm -> dense/GEMV -> GELU -> dense/GEMV -> softmax
```

The AVX2 path uses #34's non-owning spans and dot product for each dense row,
then #35's LayerNorm, activation, and stable softmax kernels. The independent
reference uses scalar double accumulation. Dimensions 5, 7, and 3 deliberately
exercise vector tails, and epsilon is explicitly `1e-5` inside LayerNorm's
square root.

## Build and run

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DDR3_BUILD_TESTS=ON -DDR3_BUILD_EXAMPLES=ON
cmake --build build --config Release --target SmallBatchInference
./build/SmallBatchInference/SmallBatchInference --self-test
./build/SmallBatchInference/SmallBatchInference
```

The self-test checks the committed fixture, scalar/AVX2 agreement for batches
1–16, softmax sums, tails, deterministic output, and rejected dimensions and
batch bounds. No Python is required.

## Local measurement

| Path | Latency (us / 2-sample batch) | Checksum |
| --- | ---: | ---: |
| Scalar double reference | 0.279 | 7.046555 |
| AVX2 `VecF8F` | 0.539 | 7.046554 |

Median of seven 200,000-iteration runs on an Intel Core i7-12700K with GCC
13.3.0, Linux x86-64, AVX2, and a Release build. This intentionally tiny,
tail-heavy fixture exposes fixed SIMD setup costs; the result is not a
universal speedup claim. Run the binary to measure the current machine.

This is an inference example, not a tensor framework or model runtime. Training,
autograd, GEMM/BLAS replacement, GPU support, and an AVX-512 path are non-goals;
AVX2 remains the demonstrated path unless the maintainers request otherwise.
