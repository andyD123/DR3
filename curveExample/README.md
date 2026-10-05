# Vector-valued curves: a QuantLib on-ramp

Keep existing QuantLib curves and `FixedRateBond` instruments. Prepare discount-factor
vectors once per distinct date, then price with the DR3 vector operations. No curve
interpolation is reimplemented by the QuantLib adapter and no analytic risk scan is
required. Each scenario may even have a different interior pillar grid.

```cpp
const auto snapshot = dr3_curve::quantlib::snapshotCurves(curves, bonds);
const auto pricer = dr3_curve::quantlib::adaptBond(*bonds.front(), snapshot);
const auto prices = pricer.prices();
// prices.npv[s]     currency NPV at the curve reference date
// prices.dirty[s]   settlement dirty price per 100 face
// prices.clean[s]   settlement clean price per 100 face
```

`curves` is a vector of ordinary QuantLib `Handle<YieldTermStructure>` objects;
`bonds` is a vector of QuantLib shared pointers to `FixedRateBond`. For a portfolio,
adapt each bond once, construct one `Prices` scratch buffer, and call `priceInto`.
The hot pricer performs no curve searches or QuantLib virtual calls. The snapshot
contains one independently allocated contiguous scenario row per date, not one flat
allocation for the entire matrix. Prepared cash flows address those rows by index.

## Build and run

The directory builds standalone or through the repository root. QuantLib is optional;
the base library does not acquire a mandatory QuantLib dependency.

```sh
# Ubuntu/Debian, for the optional QuantLib example:
sudo apt-get install libquantlib0-dev pkg-config

cmake -S curveExample -B build-curve -DCMAKE_BUILD_TYPE=Release \
  -DDR3_CURVE_ISA=AVX2 -DDR3_CURVE_ENABLE_QUANTLIB=ON
cmake --build build-curve --parallel 2
ctest --test-dir build-curve --output-on-failure
./build-curve/quantlib_curve_demo --bonds 2000 --scenarios 65 --repeats 3
./build-curve/curveExample --bonds 2000 --scenarios 65 --repeats 3
```

The ordinary root build now includes native curve regression tests whenever
`DR3_BUILD_TESTS=ON`, and the native demo whenever `DR3_BUILD_EXAMPLES=ON`.
Enable `DR3_CURVE_ENABLE_QUANTLIB=ON` at the root to include the adapter too:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DDR3_BUILD_TESTS=ON \
  -DDR3_BUILD_EXAMPLES=OFF -DDR3_CURVE_ENABLE_QUANTLIB=ON
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

A fresh root build defaults `DR3_CURVE_ISA` to `DR3_ISA`; the standalone default
remains AVX2. Disabling tests or examples suppresses those corresponding targets
rather than silently re-enabling them. QuantLib is never fetched automatically.

On multi-configuration generators use `--config Release` to build, `-C Release`
with CTest, and executables under `build-curve/Release/`. A CMake/pkg-config-visible
QuantLib installation is needed when that option is enabled. `SSE2` and `AVX512`
are alternate build targets; only execute an ISA supported by the host CPU.

Dependency-free scalar regression tests can also run from this directory alone:

```sh
cmake -S curveExample -B build-core -DDR3_CURVE_CORE_ONLY=ON
cmake --build build-core
ctest --test-dir build-core --output-on-failure
```

For new-target memory/UB checking, configure a second build with
`-DDR3_CURVE_SANITIZERS=ON -DCMAKE_BUILD_TYPE=Debug`. The scoped CI disables leak
reporting for the inherited DR3 allocation-pool lifetime policy; address and
undefined-behaviour checking remain enabled. The dependency-free core tests can
run with leak detection enabled. The scoped CI now also builds and runs the
QuantLib adapter and demo with ASan/UBSan instrumentation. The installed QuantLib
shared library itself is not rebuilt or instrumented.

## What the demonstration measures

The deterministic portfolio has varying coupons, notionals, maturities and initial
stub periods. QuantLib produces the schedules, amounts and settlement conventions.
Scenario curves include native `ZeroCurve` and `ZeroSpreadedTermStructure` objects.
For **every** bond and scenario the harness reconciles NPV, dirty price and clean
price against `DiscountingBondEngine(curve, false)`, before and after timing.
It forces recalculation rather than timing QuantLib `LazyObject` cache hits.

Three timings are reported: normal QuantLib engine repricing, scalar pricing of
prepared cached data, and DR3 pricing of that same prepared data. All write the
full NPV/dirty/clean result set. Snapshot construction and pricer preparation have
separate timings. Their sum with median pricing is labelled a sum of cold-start
components, not a direct measurement of cold-run latency. Input curve and bond
construction is shared and excluded from all three pricing timings.

This separates date/value reuse and convention extraction from vectorisation.
It is not a claim that all improvement is SIMD, or that all QuantLib workloads
receive the same improvement. The reported FNV64 hashes canonical little-endian
computed binary64 values only; it does not hash filenames, code or executables.
FMA can change the last bits relative to a scalar sum, so reconciliation uses a
stated numerical tolerance (`demo_support.h`) rather than claiming bit identity.

## Contracts and deliberate limits

- The first QuantLib on-ramp supports **fixed-rate bonds**, not arbitrary floating,
  callable or option-bearing instruments. A separate native vector test exercises
  a scenario-dependent floating coupon with distinct projection/discount curves;
  it is not a validated QuantLib floating-rate-bond adapter.
- Curves must have the same reference date, which must equal the evaluation date.
  Reference-date flows are excluded in the benchmark. Settlement-date exclusions,
  ex-coupon decisions, redemption amounts, notional and accrued interest come
  from QuantLib. `Prices` reports NPV and per-100 prices distinctly.
- A snapshot is an explicitly **frozen valuation**. It owns its discount values.
  Updating a quote, relinking a curve, rolling the date, changing an instrument,
  or changing the effective `includeTodaysCashFlows` setting requires rebuilding
  the snapshot AND re-adapting the instruments. Existing pricers intentionally retain their old valuation and have no observer hooks.
  Mixing an old snapshot with a new evaluation-date cash-flow setting during
  adaptation is rejected. QuantLib's explicit override for today's cash flows
  is honoured; the benchmark deliberately sets it to false.
- Snapshot creation and QuantLib access require a stable market/settings state.
  Snapshot reads can be shared, but concurrent DR3 allocation is subject to the
  underlying runtime contract; this example does not establish multithread scaling.
- Each source curve retains its own extrapolation policy. Missing support for a
  required maturity raises the source QuantLib exception rather than silently
  extending the curve.

## Curve-branch repairs

The scalar/vector native `Curve` and `Curve2` APIs originated in branch `Curve`
(`44e0a89a`). They are repaired and integrated against current DR3 `main`, rather
than replacing its newer vector and allocator fixes with the old branch.

Endpoint queries no longer access `pos+1` past the last pillar. Integer-date
interpolation forms overflow-safe unsigned differences before floating conversion,
so adjacent 64-bit dates above 2^53 do not collapse. Floating-date intervals with
opposite-sign endpoints whose difference overflows are scaled before division.
Invalid, empty, duplicated, unordered, nonfinite or width-mismatched inputs are rejected without replacing a valid curve.
Single-pillar curves work. Cache updates invalidate interpolation and extrapolation
results. LRU copies rebuild their iterator index. `valueAtRef` exposes cache hits
without vector copies, with explicit eviction/reset lifetimes; `valueAt` retains
the original owning-return interface.

**Numerical-definition correction:** `ZeroInterpCached` now caches precisely the
same linear-zero interpolation as `ZeroInterp`. The old undocumented `MyCalc`
forward recurrence, inconsistent initialisation, missing forward integration and
day/year ambiguity are not retained. Annualised zero rates take year fractions;
zero extrapolation is a flat boundary **rate converted to a discount factor**,
not a raw zero rate returned as a discount factor. This changes results from the
broken prototype and is covered by scalar-reference tests. The QuantLib adapter
is independent of this native zero-curve choice.

The scratch benchmark which reused one fixed date's discount vector for every
cash flow is replaced: each prepared payment now records its own correct date
row. No pointers into temporary vectors and no raw result-array allocation remain.

## Evidence

`curve_core_tests`: endpoints, validation, integer interpolation, cache lifecycle,
copy/move, zero/negative rates and extrapolation.

`curve_edge_tests`: adjacent signed/unsigned 64-bit dates, signed range-crossing,
extreme finite floating dates and subnormal intervals, with cached/uncached parity.

`curve_vector_tests`: scalar-reference reconciliation at scenario widths 1, 3, 7,
8, 9, 17, 65 and 200; fixed and floating coupon expressions; frozen rows and invalid
sizes/discounts.

`quantlib_adapter_tests`: full bond/scenario price comparisons, reference-date and
settlement-date flows, ex-coupon dates, expired bonds, different curve types/grids,
negative rates, quote refresh/frozen snapshots, today's cash-flow setting changes,
and invalid adapter inputs. Every output vector is checked before writing any
prices, including independently resized/scalar buffers; rejected buffers retain
their original contents.

A test definition is not evidence of a passing run. See the scoped CI logs and
attached run report for the commit actually executed.
