#include <stdio.h>
#include <iostream>
#include <math.h>



double normalCDF(double value)
{
    return 0.5 * erfc(-value * std::sqrt(0.5));
}

double
phi(double x)
{
    static const double RT2PI = sqrt(4.0 * acos(0.0));

    static const double SPLIT = 7.07106781186547;

    static const double N0 = 220.206867912376;
    static const double N1 = 221.213596169931;
    static const double N2 = 112.079291497871;
    static const double N3 = 33.912866078383;
    static const double N4 = 6.37396220353165;
    static const double N5 = 0.700383064443688;
    static const double N6 = 3.52624965998911e-02;
    static const double M0 = 440.413735824752;
    static const double M1 = 793.826512519948;
    static const double M2 = 637.333633378831;
    static const double M3 = 296.564248779674;
    static const double M4 = 86.7807322029461;
    static const double M5 = 16.064177579207;
    static const double M6 = 1.75566716318264;
    static const double M7 = 8.83883476483184e-02;

    const double z = fabs(x);
    double c = 0.0;

    if (z <= 37.0)
    {
        const double e = exp(-z * z / 2.0);
        if (z < SPLIT)
        {
            const double n = (((((N6 * z + N5) * z + N4) * z + N3) * z + N2) * z + N1) * z + N0;
            const double d = ((((((M7 * z + M6) * z + M5) * z + M4) * z + M3) * z + M2) * z + M1) * z + M0;
            c = e * n / d;
        }
        else
        {
            const double f = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));

            //auto  d = (z + 2.0) * (z + 4.0) == (z + 2.) * z + z * 4. + 8. = ((z + 6.) * z + 8.);
            //auto  num =  (z + 3.) * (z + 13. / 20.);   = (z +73./20.)*z +39./20)

           // auto d = (z + 2.) * (z + 4.0);//((z + 6.) * z + 8.);
           // auto n = (z + 3.) * (z + 13. / 20.);///(z + 73. / 20.) * z + 39. / 20.;

           // auto f_dash = z+ n / d;

            c = e / (RT2PI * f);
           // c = e /(f_dash*RT2PI);
        }
    }
    return x <= 0.0 ? c : 1 - c;
}



int main(int c, char** argv)
{
    for (double x = 0.01; x < 8.0; x += 0.01)
    {
        std::cout << " blah \n" << x << "\n" << phi(x) << "\n" << normalCDF(x) << "\n" << phi(x) - normalCDF(x) << "\n";
    }



	return 0;
}