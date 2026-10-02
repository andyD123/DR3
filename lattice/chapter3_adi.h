#pragma once

struct Chapter3AdiResult
{
    double center = 0.0;
    double exact_center = 0.0;
    double max_abs_error = 0.0;
};

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
