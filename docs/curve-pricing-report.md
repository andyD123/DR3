# QuantLib-to-DR3 curve/pricer on-ramp: executed results

Date: 5 October 2026.

Tested implementation commit: `859a6e8d284927aed52fa292fcee720a78be8c4d`.
Base: DR3 main `5518aba17607daa0fe1e3de8f9bdaea6fa56b045`.
Branch: `fix/curve-pricing-tools`.

[Successful scoped CI run](https://github.com/andyD123/DR3/actions/runs/37375077929)
(job `111981569780`). This report is a documentation-only follow-up to that
implementation commit. Main and the historical Curve branch were not modified.

## The on-ramp

```cpp
const auto snapshot = dr3_curve::quantlib::snapshotCurves(curves, bonds);
const auto pricer = dr3_curve::quantlib::adaptBond(*bonds.front(), snapshot);
const auto prices = pricer.prices();
```

The inputs are existing QuantLib curve handles and FixedRateBond objects.
QuantLib still evaluates discount factors and supplies cash-flow amounts,
payment dates, settlement, notional and accrued interest. Scenario curves include
ZeroCurve and ZeroSpreadedTermStructure, with different interior pillar grids.
No analytic reverse risk scan is used. See [the example README](../curveExample/README.md)
for the API, build instructions and explicit snapshot lifetime/refresh rules.

## Measured QuantLib comparison

Environment: Ubuntu 24.04.5 GitHub-hosted VM, AMD EPYC 7763 (4 virtual CPUs
visible), GCC 13.3.0, QuantLib 1.33, CMake 3.31.6, Release AVX2 / four doubles.
The pricing loops are single-threaded; this is not a thread-scaling measurement.

Command:

```sh
./build-curve/quantlib_curve_demo --bonds 2000 --scenarios 65 --repeats 3
```

The deterministic fixture contains 2,000 fixed-rate bonds, 65 scenarios and
86 distinct requested dates. It produces 130,000 bond/scenario valuations and
390,000 numeric outputs (NPV, dirty price and clean price). The strong date reuse
is an important characteristic of this synthetic workload, not a universal
property of portfolios.

| Measurement | Milliseconds |
| --- | ---: |
| QuantLib DiscountingBondEngine repricing | 2294.310419 |
| Scalar pricing of prepared cached data | 4.969744 |
| DR3 pricing of the same prepared data | 1.194923 |
| One-time QuantLib discount snapshot | 6.871669 |
| One-time instrument adaptation | 25.114804 |
| Snapshot + adaptation + median DR3 pricing | 33.181396 |

The three pricing timings are medians of three repeats. The QuantLib path
explicitly recalculates each bond, including the one-scenario case, rather than
timing LazyObject cache hits. All three paths write the full NPV/dirty/clean
result arrays. Source curve and bond construction is shared and excluded.

The DR3/scalar-cached comparison is 4.159x. The warm full-QuantLib/DR3 ratio is
1920.05x, but includes removal of repeated curve lookup, financial-convention
work and object-level pricing overhead; it is NOT a 1920x SIMD claim. The ratio
using the sum of preparation and median pricing is 69.14x. That sum is explicitly
labelled a sum of separately measured cold-start components, NOT a direct
end-to-end cold-run latency. Three repeats on one hosted VM are demonstration
evidence, not a comprehensive performance study or a universal QuantLib claim.

## Numerical reconciliation

Every result is checked against QuantLib before and after the timed loops.
Observed maximum absolute differences in the main benchmark:

| Output | Maximum absolute difference |
| --- | ---: |
| NPV in currency units | 4.547473509e-13 |
| Dirty price per 100 face | 5.684341886e-14 |
| Clean price per 100 face | 5.684341886e-14 |

Reconciliation PASS. The configured test bound is
`2e-11 * (1 + abs(reference))`, with nonfinite values rejected.

The first bond's base-scenario NPV is 992.904712593. Displayed base dirty and
clean prices are respectively 99.30135897 and 99.29036996.
Computed-output FNV64: `0xc5b2a58d3f678a21`.
The hash uses only canonical little-endian IEEE-754 binary64 result bytes, not
filenames, source files, executable bytes or metadata. FMA and differing
implementations can change low bits; cross-ISA bit identity is not claimed.

## Executed validation

| Validation | Result |
| --- | --- |
| Existing DR3 CTest regression entries | 3/3 PASS |
| Native scalar curve checks | 78 PASS |
| Native vector curve/pricing checks | 2,857 PASS |
| QuantLib adapter checks | 2,740 PASS |
| New Release CTest entries, including QuantLib | 6/6 PASS |
| New native Debug ASan/UBSan CTest entries | 4/4 PASS |

The existing regression configuration retains main's pre-existing exclusion of
TestFilterSelect.ApplyFilterB. No additional existing tests were disabled.
The native sanitizer run has leak reporting disabled for DR3's inherited pool
lifetime policy; address and undefined-behaviour checks remain active. QuantLib
adapter tests ran in Release; this CI run did not instrument the QuantLib adapter
or the installed QuantLib library with sanitizers.

Native vector coverage includes scenario widths 1, 3, 7, 8, 9, 17, 65 and 200,
independent scalar references, fixed cash flows and a scenario-dependent floating
coupon using distinct projection and discount curves. The latter is not a
QuantLib floating-rate-bond adapter.

QuantLib tests include actual reference-date and settlement-date cash flows,
ex-coupon periods, expired bonds, different curve implementations/grids, negative
rates, old frozen snapshots after quote changes, fresh-snapshot reconciliation,
and rejected invalid adapter inputs.

Additional local native-only builds on an AMD EPYC 9V74 virtual environment
passed all four native CTest entries with GCC 14.2 Release SSE2, AVX2 and AVX512,
Clang 17 Release AVX2, and GCC Debug AVX2 ASan/UBSan. These local builds do not
constitute additional QuantLib-version coverage. No Windows/MSVC, ARM or Metal
validation was performed for the new adapter in this exercise.

## Native Curve-branch repairs and limits

The old Curve branch was not merged wholesale. The curve abstraction was brought
onto current main's DR3 runtime, preserving its newer allocator/compiler fixes.
The replacement curve tool fixes final-pillar access, integer-date interpolation,
input validation, cache invalidation, LRU iterator-index copying and reference
lifetimes. Its replacement benchmark uses every cash flow's actual date row,
not the old single fixed discount vector for all dates.

There is a deliberate mathematical-definition correction: native ZeroInterpCached
now caches the same linear-zero-rate interpolation as ZeroInterp. The broken,
undocumented MyCalc forward recurrence is not preserved. Native zero curves use
annualised continuously compounded rates and year fractions; boundary zero-rate
extrapolation is converted to discount factors. The QuantLib adapter does not use
this native interpolation definition: it retains each source curve's own logic.

The QuantLib on-ramp currently supports fixed-rate bonds only. It requires common
curve reference dates equal to the evaluation date. Snapshots and adapted pricers
represent an explicitly frozen valuation. Quote changes, curve relinking, date
rolls or instrument changes require rebuilding the snapshot AND re-adapting the
instruments; existing objects intentionally retain the old valuation. QuantLib
remains an optional example dependency, not a new mandatory DR3 dependency.
