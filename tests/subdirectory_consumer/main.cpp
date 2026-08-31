#include <VecX/dr3.h>

#include <vector>

int main()
{
    std::vector<double> source{1.0, 2.0, 3.0, 4.0, 5.0};
    DRC::VecD4D::VecXX input(source);
    auto doubleValue = [](auto value) { return value + value; };
    auto output = transform(doubleValue, input);

    return output.size() == 5 && output[0] == 2.0 && output[4] == 10.0
        ? 0
        : 1;
}
