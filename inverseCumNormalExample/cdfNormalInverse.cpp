#include "cdfNormalInverse.h"


#include <algorithm>
#include <random>
#include <numeric>
#include <iterator>
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>  




//#include "../Vectorisation/VecX/curve.h"
#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/apply_operation.h"
#include "../Vectorisation/VecX/vec_d.h"
#include "../Vectorisation/VecX/vec_bool.h"
#include "../Vectorisation/VecX/vec_view.h"

#include "../Vectorisation/VecX/target_name_space.h"

#include <immintrin.h>

//using namespace DRC::VecDb;
//using namespace DRC::VecD2D;
//using namespace DRC::VecD4D;
//using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;



#include "cdfNormalInverse.h"

//#include "norm.h"






#include <chrono>

#include <iostream>
#include <functional>

double getnull(double)
{
	return 0.0;
}



//double sVML_invCumNormal(double x);



// @WichuraQuantile   
// see R implemerntation
long double qnorm8(long double p)
{

	static long double a[] = { 2509.0809287301226727 , 33430.575583588128105, 67265.770927008700853, 45921.953931549871457, 13731.693765509461125,  1971.5909503065514427, 133.14166789178437745,3.387132872796366608 };
	static long double b[] = { 5226.495278852854561, 28729.085735721942674,    39307.89580009271061, 21213.794301586595867, 5394.1960214247511077,   687.1870074920579083, 42.313330701600911252 };
	static long double c[] = { 7.7454501427834140764e-4 , .0227238449892691845833 ,.24178072517745061177, 1.27045825245236838258 ,  3.64784832476320460504, 5.7694972214606914055, 4.6303378461565452959, 1.42343711074968357734 };
	static long double d[] = { 1.05075007164441684324e-9 , 5.475938084995344946e-4, .0151986665636164571966, .14810397642748007459, .68976733498510000455,  1.6763848301838038494,  2.05319162663775882187,1. };
	static long double e[] = { 2.01033439929228813265e-7 ,   2.71155556874348757815e-5,   .0012426609473880784386, .026532189526576123093, .29656057182850489123,   1.7848265399172913358, 5.4637849111641143699, 6.6579046435011037772 };
	static long double f[] = { 2.04426310338993978564e-15 , 1.4215117583164458887e-7, 1.8463183175100546818e-5,  7.868691311456132591e-4, .0148753612908506148525,.13692988092273580531, .59983220655588793769, 1. };


	long double val = 0.0;
	long double   q = p - 0.5;


	if (fabs(q) <= .425)
	{
		long double r = .180625 - q * q;
		val =

			q * (((((((r * a[0] +
				a[1]) * r + a[2]) * r +
				a[3]) * r + a[4]) * r +
				a[5]) * r + a[6]) * r +
				a[7])
			/ (((((((r * b[0] +
				b[1]) * r + b[2]) * r +
				b[3]) * r + b[4]) * r +
				b[5]) * r + b[6]) * r + 1.);

	}
	else
	{ // closer than 0.075 from {0,1} boundary 

	   // r = min(p, 1-p) < 0.075 
		long double  r = std::min(p, 1 - p);
		r = sqrt(-log(r));



		// <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 
		if (r <= 5.)
		{
			r += -1.6;
			val = (((((((r * c[0] + c[1]) * r + c[2]) * r + c[3]) * r + c[4]) * r + c[5]) * r + c[6]) * r + c[7])
				/ (((((((r * d[0] + d[1]) * r + d[2]) * r + d[3]) * r + d[4]) * r + d[5]) * r + d[6]) * r + 1.);


		}
		else
		{ // very close to  0 or 1 
			r += -5.;

			val = (((((((r * e[0] + e[1]) * r + e[2]) * r + e[3]) * r + e[4]) * r + e[5]) * r + e[6]) * r + e[7])
				/ (((((((r * f[0] + f[1]) * r + f[2]) * r + f[3]) * r + f[4]) * r + f[5]) * r + f[6]) * r + 1.);

		}

		long double valMult = (q < 0.0) ? -1.0 : 1.0;
		val *= valMult;
	}

	return val;
}




