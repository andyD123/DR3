#include "chapter3_examples.h"
#include "chapter3_adi.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <vector>

namespace
{
constexpr double kS = 100.0;
constexpr double kK = 100.0;
constexpr double kT = 1.0;
constexpr double kSigma = 0.2;
constexpr double kRate = 0.06;
constexpr double kDiv = 0.03;
constexpr int kN = 3;
constexpr int kNj = 3;
constexpr double kDx = 0.2;

std::vector<double> chapter3AssetGrid()
{
    std::vector<double> stock(2 * kNj + 1);
    for (int j = -kNj; j <= kNj; ++j)
        stock[j + kNj] = kS * std::exp(j * kDx);
    return stock;
}

double chapter3TrinomialEuropeanCall()
{
    const double dt = kT / kN;
    const double nu = kRate - kDiv - 0.5 * kSigma * kSigma;

    const double pu = 0.5 * (
        (dt * kSigma * kSigma + nu * nu * dt * dt) / (kDx * kDx)
        + nu * dt / kDx);
    const double pd = 0.5 * (
        (dt * kSigma * kSigma + nu * nu * dt * dt) / (kDx * kDx)
        - nu * dt / kDx);
    const double pm = 1.0 -
        (dt * kSigma * kSigma + nu * nu * dt * dt) / (kDx * kDx);
    const double disc = std::exp(-kRate * dt);

    const auto stock = chapter3AssetGrid();
    std::vector<double> next(stock.size());
    std::transform(stock.begin(), stock.end(), next.begin(),
        [](double s) { return std::max(0.0, s - kK); });

    for (int i = kN - 1; i >= 0; --i)
    {
        auto current = next;
        for (int j = -i; j <= i; ++j)
        {
            const int p = j + kNj;
            current[p] = disc * (
                pu * next[p + 1] +
                pm * next[p] +
                pd * next[p - 1]);
        }
        next.swap(current);
    }

    return next[kNj];
}

struct ExplicitCoefficients
{
    double pu;
    double pm;
    double pd;
};

ExplicitCoefficients explicitCoefficients()
{
    const double dt = kT / kN;
    const double nu = kRate - kDiv - 0.5 * kSigma * kSigma;
    return {
        0.5 * dt * (kSigma * kSigma / (kDx * kDx) + nu / kDx),
        1.0 - dt * kSigma * kSigma / (kDx * kDx) - kRate * dt,
        0.5 * dt * (kSigma * kSigma / (kDx * kDx) - nu / kDx)
    };
}

double chapter3ExplicitEuropeanCall()
{
    const auto c = explicitCoefficients();
    const auto stock = chapter3AssetGrid();

    std::vector<double> next(stock.size());
    std::transform(stock.begin(), stock.end(), next.begin(),
        [](double s) { return std::max(0.0, s - kK); });

    const int J = 2 * kNj;
    for (int i = kN - 1; i >= 0; --i)
    {
        std::vector<double> current(next.size(), 0.0);
        for (int j = 1; j < J; ++j)
            current[j] = c.pu * next[j + 1] + c.pm * next[j] + c.pd * next[j - 1];

        current[0] = current[1];
        current[J] = current[J - 1] + stock[J] - stock[J - 1];
        next.swap(current);
    }

    return next[kNj];
}

double chapter3ExplicitAmericanPut()
{
    const auto c = explicitCoefficients();
    const auto stock = chapter3AssetGrid();

    std::vector<double> exercise(stock.size());
    std::transform(stock.begin(), stock.end(), exercise.begin(),
        [](double s) { return std::max(0.0, kK - s); });

    auto next = exercise;
    const int J = 2 * kNj;

    for (int i = kN - 1; i >= 0; --i)
    {
        std::vector<double> current(next.size(), 0.0);
        for (int j = 1; j < J; ++j)
            current[j] = c.pu * next[j + 1] + c.pm * next[j] + c.pd * next[j - 1];

        current[0] = current[1] + stock[1] - stock[0];
        current[J] = current[J - 1];

        for (int j = 0; j <= J; ++j)
            current[j] = std::max(current[j], exercise[j]);

        next.swap(current);
    }

    return next[kNj];
}

struct ImplicitCoefficients
{
    double pu;
    double pm;
    double pd;
};

ImplicitCoefficients implicitCoefficients(double theta)
{
    const double dt = kT / kN;
    const double nu = kRate - kDiv - 0.5 * kSigma * kSigma;

    return {
        -theta * dt * (kSigma * kSigma / (kDx * kDx) + nu / kDx),
        1.0 + 2.0 * theta * dt * kSigma * kSigma / (kDx * kDx)
            + 2.0 * theta * kRate * dt,
        -theta * dt * (kSigma * kSigma / (kDx * kDx) - nu / kDx)
    };
}

void solveImplicitStep(
    const std::vector<double>& rhs,
    std::vector<double>& result,
    double pu,
    double pm,
    double pd,
    double lambdaL,
    double lambdaU)
{
    const int J = static_cast<int>(rhs.size()) - 1;
    std::vector<double> pmp(rhs.size(), 0.0);
    std::vector<double> pp(rhs.size(), 0.0);

    pmp[1] = pm + pd;
    pp[1] = rhs[1] + pd * lambdaL;

    for (int j = 2; j < J; ++j)
    {
        pmp[j] = pm - pu * pd / pmp[j - 1];
        pp[j] = rhs[j] - pp[j - 1] * pd / pmp[j - 1];
    }

    result[J] = (pp[J - 1] + pmp[J - 1] * lambdaU) /
        (pu + pmp[J - 1]);
    result[J - 1] = result[J] - lambdaU;

    for (int j = J - 2; j != 0; --j)
        result[j] = (pp[j] - pu * result[j + 1]) / pmp[j];

    result[0] = result[1] - lambdaL;
}

double chapter3ImplicitAmericanPut()
{
    const auto c = implicitCoefficients(0.5);
    const auto stock = chapter3AssetGrid();

    std::vector<double> exercise(stock.size());
    std::transform(stock.begin(), stock.end(), exercise.begin(),
        [](double s) { return std::max(0.0, kK - s); });

    auto next = exercise;
    const double lambdaL = -(stock[1] - stock[0]);
    const double lambdaU = 0.0;

    for (int i = kN - 1; i >= 0; --i)
    {
        std::vector<double> current(next.size(), 0.0);
        solveImplicitStep(next, current, c.pu, c.pm, c.pd, lambdaL, lambdaU);

        for (std::size_t j = 0; j < current.size(); ++j)
            current[j] = std::max(current[j], exercise[j]);

        next.swap(current);
    }

    return next[kNj];
}

double chapter3CrankNicolsonAmericanPut()
{
    const auto c = implicitCoefficients(0.25);
    const auto stock = chapter3AssetGrid();

    std::vector<double> exercise(stock.size());
    std::transform(stock.begin(), stock.end(), exercise.begin(),
        [](double s) { return std::max(0.0, kK - s); });

    auto next = exercise;
    const int J = 2 * kNj;
    const double lambdaL = -(stock[1] - stock[0]);
    const double lambdaU = 0.0;

    for (int i = kN - 1; i >= 0; --i)
    {
        std::vector<double> transformed(next.size(), 0.0);
        transformed[1] =
            -c.pu * next[2] -
            (c.pm - 2.0) * next[1] -
            c.pd * next[0] +
            c.pd * lambdaL;

        for (int j = 2; j < J; ++j)
            transformed[j] =
                -c.pu * next[j + 1] -
                (c.pm - 2.0) * next[j] -
                c.pd * next[j - 1];

        // Reuse the implicit elimination with the Crank-Nicolson right-hand side.
        std::vector<double> pmp(next.size(), 0.0);
        std::vector<double> pp(next.size(), 0.0);
        pmp[1] = c.pm + c.pd;
        pp[1] = transformed[1];

        for (int j = 2; j < J; ++j)
        {
            pmp[j] = c.pm - c.pu * c.pd / pmp[j - 1];
            pp[j] = transformed[j] - pp[j - 1] * c.pd / pmp[j - 1];
        }

        std::vector<double> current(next.size(), 0.0);
        current[J] = (pp[J - 1] + pmp[J - 1] * lambdaU) /
            (c.pu + pmp[J - 1]);
        current[J - 1] = current[J] - lambdaU;

        for (int j = J - 2; j != 0; --j)
            current[j] = (pp[j] - c.pu * current[j + 1]) / pmp[j];

        current[0] = current[1] - lambdaL;

        for (int j = 0; j <= J; ++j)
            current[j] = std::max(current[j], exercise[j]);

        next.swap(current);
    }

    return next[kNj];
}

bool nearPrinted(double value, double printed)
{
    // Chapter 3 reports the worked examples to four decimal places.
    return std::abs(value - printed) <= 5.0e-4;
}

void printCheck(
    std::ostream& out,
    const char* name,
    double value,
    double expected,
    bool pass)
{
    out << "  " << std::left << std::setw(34) << name
        << std::right << std::fixed << std::setprecision(8)
        << value << "  expected " << std::setprecision(4) << expected
        << "  " << (pass ? "PASS" : "FAIL") << "\n";
}
} // namespace

Chapter3RegressionSummary runChapter3Regression(std::ostream& out)
{
    Chapter3RegressionSummary r;

    r.trinomial_european_call = chapter3TrinomialEuropeanCall();
    r.explicit_european_call = chapter3ExplicitEuropeanCall();
    r.explicit_american_put = chapter3ExplicitAmericanPut();
    r.implicit_american_put = chapter3ImplicitAmericanPut();
    r.crank_nicolson_american_put = chapter3CrankNicolsonAmericanPut();

    const auto adiReference = chapter3AdiReference();
    const auto adiDr3 = chapter3AdiDr3();

    r.adi_reference_center = adiReference.center;
    r.adi_dr3_center = adiDr3.center;
    r.adi_exact_center = adiReference.exact_center;
    r.adi_reference_max_error = adiReference.max_abs_error;
    r.adi_dr3_max_error = adiDr3.max_abs_error;

    const bool p1 = nearPrinted(r.trinomial_european_call, 8.4253);
    const bool p2 = nearPrinted(r.explicit_european_call, 8.5455);
    const bool p3 = nearPrinted(r.explicit_american_put, 6.0058);
    const bool p4 = nearPrinted(r.implicit_american_put, 4.9221);
    const bool p5 = nearPrinted(r.crank_nicolson_american_put, 5.4184);

    const bool adiAgreement =
        std::abs(r.adi_reference_center - r.adi_dr3_center) <= 1.0e-12;
    const bool adiAccuracy = r.adi_reference_max_error <= 3.0e-4
        && r.adi_dr3_max_error <= 3.0e-4;

    out << "\nChapter 3 source-fidelity regressions\n";
    printCheck(out, "trinomial European call", r.trinomial_european_call, 8.4253, p1);
    printCheck(out, "explicit FD European call", r.explicit_european_call, 8.5455, p2);
    printCheck(out, "explicit FD American put", r.explicit_american_put, 6.0058, p3);
    printCheck(out, "implicit FD American put", r.implicit_american_put, 4.9221, p4);
    printCheck(out, "Crank-Nicolson American put", r.crank_nicolson_american_put, 5.4184, p5);

    out << "  ADI equation (3.71) manufactured mode\n"
        << "    scalar center = " << std::setprecision(12) << r.adi_reference_center << "\n"
        << "    MDSpan center = " << r.adi_dr3_center << "\n"
        << "    exact center  = " << r.adi_exact_center << "\n"
        << "    scalar max error = " << r.adi_reference_max_error << "\n"
        << "    MDSpan max error = " << r.adi_dr3_max_error << "\n"
        << "    reference/MDSpan agreement " << (adiAgreement ? "PASS" : "FAIL") << "\n"
        << "    manufactured-solution error " << (adiAccuracy ? "PASS" : "FAIL") << "\n";

    r.passed = p1 && p2 && p3 && p4 && p5 && adiAgreement && adiAccuracy;
    out << "Chapter 3 regressions: " << (r.passed ? "PASS" : "FAIL") << "\n\n";
    return r;
}
