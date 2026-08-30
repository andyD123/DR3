#include "tree_pricers.h"

#include <chrono>
#include <iomanip>
#include <iostream>

int main()
{
    using namespace dr3::lattice;

    const VanillaOption option{100.0, 100.0, 0.20, 0.05, 0.0, 1.0, OptionType::Call};
    const TreeConfig config{1024};
    const auto start = std::chrono::steady_clock::now();
    const double price = europeanTrinomial(option, config);
    const auto elapsed = std::chrono::steady_clock::now() - start;

    std::cout << std::fixed << std::setprecision(8)
              << "European call: " << price << '\n'
              << "Steps: " << config.steps << '\n'
              << "Elapsed (us): "
              << std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count()
              << '\n';
}
