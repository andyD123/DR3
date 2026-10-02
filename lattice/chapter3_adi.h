#pragma once

struct Chapter3AdiResult
{
    double center = 0.0;
    double exact_center = 0.0;
    double max_abs_error = 0.0;
};

struct Chapter3AdiTransform
{
    double lambda1 = 0.0;
    double lambda2 = 0.0;
    double e11 = 0.0;
    double e12 = 0.0;
    double e21 = 0.0;
    double e22 = 0.0;
    double alpha1 = 0.0;
    double alpha2 = 0.0;
    double a1 = 0.0;
    double a2 = 0.0;
    double a3 = 0.0;
    double y2_scale = 0.0;
    double rotated_off_diagonal = 0.0;
};

Chapter3AdiTransform chapter3AdiTransform(
    double sigma1,
    double sigma2,
    double rho,
    double nu1,
    double nu2,
    double rate);

Chapter3AdiResult chapter3AdiReference(
    int grid_points = 33,
    int time_steps = 80,
    double maturity = 0.25,
    double lambda1 = 0.2);

Chapter3AdiResult chapter3AdiDr3(
    int grid_points = 33,
    int time_steps = 80,
    double maturity = 0.25,
    double lambda1 = 0.2);
