#pragma once

#include <iosfwd>

struct Chapter3RegressionSummary
{
    bool passed = false;

    double trinomial_european_call = 0.0;
    double explicit_european_call = 0.0;
    double explicit_american_put = 0.0;
    double implicit_american_put = 0.0;
    double crank_nicolson_american_put = 0.0;

    double adi_reference_center = 0.0;
    double adi_dr3_center = 0.0;
    double adi_exact_center = 0.0;
    double adi_reference_max_error = 0.0;
    double adi_dr3_max_error = 0.0;
    double adi_transform_off_diagonal = 0.0;
};

Chapter3RegressionSummary runChapter3Regression(std::ostream& out);
