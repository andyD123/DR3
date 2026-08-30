# Vectorised Monte Carlo European call

Deterministic uniforms, inverse-normal via DR3 `cdfnorminv`, antithetic
variates, discounted call payoff, standard error, and a 95% confidence
interval. Scalar and AVX2 implementations are compared with analytic
Black-Scholes.

```
cmake -S . -B build -DDR3_BUILD_TESTS=ON -DDR3_BUILD_EXAMPLES=ON
cmake --build build --target monteCarloExample
./build/monteCarloExample/monteCarloExample --self-test
./build/monteCarloExample/monteCarloExample --paths 65536
```

Printed errors are local measurements, not a claimed speedup.
