# Vectorised portfolio stress and P&L aggregation

Prices a deterministic option book, applies spot and volatility shocks,
filters losses beyond a threshold with DR3 views, and compares naive,
pairwise, and Kahan aggregation.

Demonstration only.

```
cmake -S . -B build -DDR3_BUILD_TESTS=ON -DDR3_BUILD_EXAMPLES=ON
cmake --build build --target portfolioStressExample
./build/portfolioStressExample/portfolioStressExample --self-test
./build/portfolioStressExample/portfolioStressExample --size 1024
```

Printed residuals are local measurements, not a claimed speedup.
