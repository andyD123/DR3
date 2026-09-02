# Exact embedding similarity and top-k search

Educational example: exact (not ANN) retrieval over a contiguous float32
embedding matrix using DR3 AVX2 (`VecF8F`) transform/reduce kernels.

APIs in `Vectorisation/VecX/similarity.h`:

- `dot_product`, `squared_l2_distance`, `cosine_similarity`
- `normalize_l2_inplace`
- `top_k_inner_product`, `top_k_cosine_similarity`, `top_k_l2`

Mismatched lengths throw `std::invalid_argument`. Cosine of a zero-norm
vector is the sentinel `0` (never NaN). Duplicate scores tie-break by the
smaller original row index. Optional score thresholds drop results below
the cutoff (or above it for L2). Top-k uses a bounded heap; it does not
fully sort the corpus.

This is not a Faiss / oneDNN / ONNX replacement.

## Build and test

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DDR3_BUILD_TESTS=ON -DDR3_BUILD_EXAMPLES=ON
cmake --build build --config Release --target VectorTest EmbeddingSearch
ctest --test-dir build -C Release -R 'Similarity|Embedding' --output-on-failure
./build/EmbeddingSearch/EmbeddingSearch --self-test
```

## Benchmark

```
./build/EmbeddingSearch/EmbeddingSearch --size 4096 --dim 128 --k 10
```

Printed times depend on processor, compiler, instruction set, and build
type. The README does not claim a speedup.
