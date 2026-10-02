#include "chapter3_adi.h"

#include "../Vectorisation/VecX/dr3.h"
#include "../Vectorisation/VecX/span.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{
constexpr double kPi = 3.141592653589793238462643383279502884;

void validateAdiArguments(int gridPoints, int timeSteps, double maturity, double lambda1)
{
    if (gridPoints < 5 || (gridPoints % 2) == 0)
        throw std::invalid_argument("Chapter 3 ADI example requires an odd grid_points >= 5");
    if (timeSteps <= 0)
        throw std::invalid_argument("Chapter 3 ADI example requires time_steps > 0");
    if (maturity <= 0.0 || lambda1 <= 0.0)
        throw std::invalid_argument("Chapter 3 ADI example requires positive maturity and lambda1");
}

void solveConstantTridiagonal(
    double a,
    double b,
    double c,
    const std::vector<double>& rhs,
    std::vector<double>& solution)
{
    const int n = static_cast<int>(rhs.size());
    solution.resize(rhs.size());

    if (n == 1)
    {
        solution[0] = rhs[0] / b;
        return;
    }

    std::vector<double> cPrime(n - 1);
    std::vector<double> dPrime(n);

    cPrime[0] = c / b;
    dPrime[0] = rhs[0] / b;

    for (int i = 1; i < n; ++i)
    {
        const double denominator = b - a * cPrime[i - 1];
        if (i < n - 1)
            cPrime[i] = c / denominator;
        dPrime[i] = (rhs[i] - a * dPrime[i - 1]) / denominator;
    }

    solution[n - 1] = dPrime[n - 1];
    for (int i = n - 2; i >= 0; --i)
        solution[i] = dPrime[i] - cPrime[i] * solution[i + 1];
}

template<class Get, class Set>
void initialiseMode(int n, Get&& coordinate, Set&& set)
{
    for (int j = 0; j < n; ++j)
    {
        const double y1 = coordinate(j);
        for (int k = 0; k < n; ++k)
        {
            const double y2 = coordinate(k);
            set(j, k, std::sin(kPi * y1) * std::sin(kPi * y2));
        }
    }
}

template<class Get>
Chapter3AdiResult measureAdiResult(
    int n,
    double maturity,
    double lambda1,
    Get&& get)
{
    const double amplitude = std::exp(-lambda1 * kPi * kPi * maturity);
    const int mid = n / 2;

    double maxError = 0.0;
    for (int j = 0; j < n; ++j)
    {
        const double y1 = static_cast<double>(j) / static_cast<double>(n - 1);
        for (int k = 0; k < n; ++k)
        {
            const double y2 = static_cast<double>(k) / static_cast<double>(n - 1);
            const double exact =
                amplitude * std::sin(kPi * y1) * std::sin(kPi * y2);
            maxError = std::max(maxError, std::abs(get(j, k) - exact));
        }
    }

    return { get(mid, mid), amplitude, maxError };
}
} // namespace

Chapter3AdiTransform chapter3AdiTransform(
    double sigma1,
    double sigma2,
    double rho,
    double nu1,
    double nu2,
    double rate)
{
    if (sigma1 <= 0.0 || sigma2 <= 0.0)
        throw std::invalid_argument("Chapter 3 ADI transform requires positive volatilities");
    if (rho <= -1.0 || rho >= 1.0)
        throw std::invalid_argument("Chapter 3 ADI transform requires -1 < rho < 1");

    // Equation (3.64): symmetric covariance matrix in log-asset coordinates.
    const double a = sigma1 * sigma1;
    const double b = rho * sigma1 * sigma2;
    const double d = sigma2 * sigma2;

    const double mean = 0.5 * (a + d);
    const double radius = std::sqrt(0.25 * (a - d) * (a - d) + b * b);
    const double lambda1 = mean + radius;
    const double lambda2 = mean - radius;

    if (lambda2 <= 0.0)
        throw std::invalid_argument("Chapter 3 ADI covariance matrix is not positive definite");

    // Normalised eigenvector for lambda1; the orthogonal complement is lambda2.
    double e11 = b;
    double e12 = lambda1 - a;
    if (std::abs(e11) + std::abs(e12) < 1.0e-15)
    {
        e11 = 1.0;
        e12 = 0.0;
    }
    const double norm = std::sqrt(e11 * e11 + e12 * e12);
    e11 /= norm;
    e12 /= norm;

    const double e21 = -e12;
    const double e22 = e11;

    // Equation (3.65): drift in the uncorrelated coordinates.
    const double alpha1 = e11 * nu1 + e12 * nu2;
    const double alpha2 = e21 * nu1 + e22 * nu2;

    // Equations (3.67)-(3.70).
    const double a1 = -alpha1 / lambda1;
    const double a2 = -alpha2 / lambda2;
    const double a3 = alpha1 * alpha1 / (2.0 * lambda1)
        + alpha2 * alpha2 / (2.0 * lambda2)
        + rate;
    const double y2Scale = std::sqrt(lambda1 / lambda2);

    // E * covariance * E^T should be diagonal.
    const double rotatedOffDiagonal =
        e11 * (a * e21 + b * e22)
        + e12 * (b * e21 + d * e22);

    return {
        lambda1, lambda2,
        e11, e12, e21, e22,
        alpha1, alpha2,
        a1, a2, a3,
        y2Scale,
        rotatedOffDiagonal
    };
}

Chapter3AdiResult chapter3AdiReference(
    int gridPoints,
    int timeSteps,
    double maturity,
    double lambda1)
{
    validateAdiArguments(gridPoints, timeSteps, maturity, lambda1);

    const int n = gridPoints;
    const double dy = 1.0 / static_cast<double>(n - 1);
    const double dt = maturity / static_cast<double>(timeSteps);
    const double invDy2 = 1.0 / (dy * dy);

    // Equations (3.73)-(3.76), multiplied by -1 to give a positive diagonal.
    const double a = -lambda1 * 0.5 * invDy2;
    const double b = 2.0 / dt + lambda1 * invDy2;
    const double c = a;

    auto index = [n](int j, int k) { return static_cast<std::size_t>(j * n + k); };

    std::vector<double> u(static_cast<std::size_t>(n * n), 0.0);
    std::vector<double> half(u.size(), 0.0);
    std::vector<double> next(u.size(), 0.0);

    initialiseMode(
        n,
        [n](int p) { return static_cast<double>(p) / static_cast<double>(n - 1); },
        [&](int j, int k, double value) { u[index(j, k)] = value; });

    std::vector<double> rhs(n - 2);
    std::vector<double> solution;

    for (int step = 0; step < timeSteps; ++step)
    {
        std::fill(half.begin(), half.end(), 0.0);

        // First half-step: y1 implicit, y2' explicit, equation (3.72).
        for (int k = 1; k < n - 1; ++k)
        {
            for (int j = 1; j < n - 1; ++j)
            {
                const double d2y2 =
                    (u[index(j, k + 1)] - 2.0 * u[index(j, k)] + u[index(j, k - 1)])
                    * invDy2;
                rhs[j - 1] = 2.0 / dt * u[index(j, k)] + 0.5 * lambda1 * d2y2;
            }

            solveConstantTridiagonal(a, b, c, rhs, solution);
            for (int j = 1; j < n - 1; ++j)
                half[index(j, k)] = solution[j - 1];
        }

        std::fill(next.begin(), next.end(), 0.0);

        // Second half-step: y2' implicit, y1 explicit, equation (3.77).
        for (int j = 1; j < n - 1; ++j)
        {
            for (int k = 1; k < n - 1; ++k)
            {
                const double d2y1 =
                    (half[index(j + 1, k)] - 2.0 * half[index(j, k)] + half[index(j - 1, k)])
                    * invDy2;
                rhs[k - 1] = 2.0 / dt * half[index(j, k)] + 0.5 * lambda1 * d2y1;
            }

            solveConstantTridiagonal(a, b, c, rhs, solution);
            for (int k = 1; k < n - 1; ++k)
                next[index(j, k)] = solution[k - 1];
        }

        u.swap(next);
    }

    return measureAdiResult(
        n, maturity, lambda1,
        [&](int j, int k) { return u[index(j, k)]; });
}

Chapter3AdiResult chapter3AdiDr3(
    int gridPoints,
    int timeSteps,
    double maturity,
    double lambda1)
{
    validateAdiArguments(gridPoints, timeSteps, maturity, lambda1);

    using Scalar = typename InstructionTraits<VecXX::INS>::FloatType;
    constexpr std::size_t width = InstructionTraits<VecXX::INS>::width;
    using Layout = Layout2D<Scalar, width, 0>;

    const int n = gridPoints;
    const std::size_t padded =
        ((static_cast<std::size_t>(n) + width - 1) / width) * width;
    const std::size_t storageSize = padded * static_cast<std::size_t>(n);

    std::vector<Scalar> uStorage(storageSize, Scalar(0));
    std::vector<Scalar> halfStorage(storageSize, Scalar(0));
    std::vector<Scalar> nextStorage(storageSize, Scalar(0));

    MDSpan<Scalar, Layout> u(uStorage.data(), n, n);
    MDSpan<Scalar, Layout> half(halfStorage.data(), n, n);
    MDSpan<Scalar, Layout> next(nextStorage.data(), n, n);

    const double dy = 1.0 / static_cast<double>(n - 1);
    const double dt = maturity / static_cast<double>(timeSteps);
    const double invDy2 = 1.0 / (dy * dy);

    const double a = -lambda1 * 0.5 * invDy2;
    const double b = 2.0 / dt + lambda1 * invDy2;
    const double c = a;

    initialiseMode(
        n,
        [n](int p) { return static_cast<double>(p) / static_cast<double>(n - 1); },
        [&](int j, int k, double value) { u(j, k) = static_cast<Scalar>(value); });

    std::vector<double> rhs(n - 2);
    std::vector<double> solution;

    for (int step = 0; step < timeSteps; ++step)
    {
        std::fill(halfStorage.begin(), halfStorage.end(), Scalar(0));

        // Rows/columns are deliberately accessed through MDSpan here. This is the
        // Chapter 3 two-dimensional layout realization; later SIMD work can replace
        // the scalar line solves without changing the numerical contract.
        for (int k = 1; k < n - 1; ++k)
        {
            for (int j = 1; j < n - 1; ++j)
            {
                const double d2y2 =
                    (static_cast<double>(u(j, k + 1))
                     - 2.0 * static_cast<double>(u(j, k))
                     + static_cast<double>(u(j, k - 1))) * invDy2;
                rhs[j - 1] =
                    2.0 / dt * static_cast<double>(u(j, k))
                    + 0.5 * lambda1 * d2y2;
            }

            solveConstantTridiagonal(a, b, c, rhs, solution);
            for (int j = 1; j < n - 1; ++j)
                half(j, k) = static_cast<Scalar>(solution[j - 1]);
        }

        std::fill(nextStorage.begin(), nextStorage.end(), Scalar(0));

        for (int j = 1; j < n - 1; ++j)
        {
            for (int k = 1; k < n - 1; ++k)
            {
                const double d2y1 =
                    (static_cast<double>(half(j + 1, k))
                     - 2.0 * static_cast<double>(half(j, k))
                     + static_cast<double>(half(j - 1, k))) * invDy2;
                rhs[k - 1] =
                    2.0 / dt * static_cast<double>(half(j, k))
                    + 0.5 * lambda1 * d2y1;
            }

            solveConstantTridiagonal(a, b, c, rhs, solution);
            for (int k = 1; k < n - 1; ++k)
                next(j, k) = static_cast<Scalar>(solution[k - 1]);
        }

        // Keep the MDSpan views stable. The source-fidelity example favours
        // auditability over avoiding this copy; the optimized realization can replace it.
        std::copy(nextStorage.begin(), nextStorage.end(), uStorage.begin());
    }

    return measureAdiResult(
        n, maturity, lambda1,
        [&](int j, int k) { return static_cast<double>(u(j, k)); });
}
