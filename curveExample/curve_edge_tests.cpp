#include "curve.h"
#include <cstdint>
#include <iostream>
#include <limits>

namespace {
int checks = 0;
void near(double actual, double expected) {
    ++checks;
    if (!std::isfinite(actual) || std::abs(actual - expected) > 1e-14)
        throw std::runtime_error("Interpolation edge regression");
}
template<class T> void interval(T lo, T middle, T hi) {
    const std::vector<T> dates{lo, hi};
    const std::vector<double> values{0.0, 2.0};
    Curve<T, double> curve;
    Curve2<T, double> cached(2);
    curve.setValues(dates.begin(), dates.end(), values.begin(), values.end());
    cached.setValues(dates.begin(), dates.end(), values.begin(), values.end());
    for (int repeat = 0; repeat < 2; ++repeat) {
        near(curve.valueAt(lo), 0.0);
        near(curve.valueAt(middle), 1.0);
        near(curve.valueAt(hi), 2.0);
        near(cached.valueAt(middle), 1.0);
    }
}
}

int main() {
    try {
        using I = std::int64_t;
        using U = std::uint64_t;
        const I imin = std::numeric_limits<I>::min();
        const I imax = std::numeric_limits<I>::max();
        const U umax = std::numeric_limits<U>::max();
        // Adjacent large integers must not collapse to the same floating value.
        interval<I>(imax - 2, imax - 1, imax);
        interval<I>(imin, imin + 1, imin + 2);
        interval<U>(umax - 2, umax - 1, umax);
        interval<I>((I(1) << 53) + 1, (I(1) << 53) + 2, (I(1) << 53) + 3);
        // Opposite-sign endpoints require an unsigned difference, not signed UB.
        interval<I>(imin, 0, imax);
        interval<std::int32_t>(-std::numeric_limits<std::int32_t>::max(), 0,
                               std::numeric_limits<std::int32_t>::max());
        const double huge = std::numeric_limits<double>::max();
        interval<double>(-huge, 0.0, huge);
        const double tiny = std::numeric_limits<double>::denorm_min();
        interval<double>(0.0, tiny, tiny * 2.0);
        std::cout << "curve_edge_tests PASS checks=" << checks << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL after " << checks << " checks: " << e.what() << '\n';
        return 1;
    }
}
