/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998       Ross Ihaka
 *  Copyright (C) 2000--2005 The R Core Team
 *  based on AS 111 (C) 1977 Royal Statistical Society
 *  and   on AS 241 (C) 1988 Royal Statistical Society
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *	double qnorm5(double p, double mu, double sigma,
 *		      int lower_tail, int log_p)
 *            {qnorm (..) is synonymous and preferred inside R}
 *
 *  DESCRIPTION
 *
 *	Compute the quantile function for the normal distribution.
 *
 *	For small to moderate probabilities, algorithm referenced
 *	below is used to obtain an initial approximation which is
 *	polished with a final Newton step.
 *
 *	For very large arguments, an algorithm of Wichura is used.
 *
 *  REFERENCE
 *
 *	Beasley, J. D. and S. G. Springer (1977).
 *	Algorithm AS 111: The percentage points of the normal distribution,
 *	Applied Statistics, 26, 118-121.
 *
 *      Wichura, M.J. (1988).
 *      Algorithm AS 241: The Percentage Points of the Normal Distribution.
 *      Applied Statistics, 37, 477-484.
 */

//#include "nmath.h"
//#include "dpq.h"

//////////////////////////////////////////
#include <math.h>
#include <algorithm>

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2000--2014 The  R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */
 /* Utilities for `dpq' handling (density/probability/quantile) */

/* give_log in "d";  log_p in "p" & "q" : */
#define give_log log_p
                            /* "DEFAULT" */
                            /* --------- */
#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)		/* 1 */
#define R_D_half (log_p ? -M_LN2 : 0.5)		// 1/2 (lower- or upper tail)


/* Use 0.5 - p + 0.5 to perhaps gain 1 bit of accuracy */
#define R_D_Lval(p)	(lower_tail ? (p) : (0.5 - (p) + 0.5))	/*  p  */
#define R_D_Cval(p)	(lower_tail ? (0.5 - (p) + 0.5) : (p))	/*  1 - p */

#define R_D_val(x)	(log_p	? log(x) : (x))		/*  x  in pF(x,..) */
#define R_D_qIv(p)	(log_p	? exp(p) : (p))		/*  p  in qF(p,..) */
#define R_D_exp(x)	(log_p	?  (x)	 : exp(x))	/* exp(x) */
#define R_D_log(p)	(log_p	?  (p)	 : log(p))	/* log(p) */
#define R_D_Clog(p)	(log_p	? log1p(-(p)) : (0.5 - (p) + 0.5)) /* [log](1-p) */

// log(1 - exp(x))  in more stable form than log1p(- R_D_qIv(x)) :
#define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))

/* log(1-exp(x)):  R_D_LExp(x) == (log1p(- R_D_qIv(x))) but even more stable:*/
#define R_D_LExp(x)     (log_p ? R_Log1_Exp(x) : log1p(-x))

#define R_DT_val(x)	(lower_tail ? R_D_val(x)  : R_D_Clog(x))
#define R_DT_Cval(x)	(lower_tail ? R_D_Clog(x) : R_D_val(x))

/*#define R_DT_qIv(p)	R_D_Lval(R_D_qIv(p))		 *  p  in qF ! */
#define R_DT_qIv(p)	(log_p ? (lower_tail ? exp(p) : - expm1(p)) \
			       : R_D_Lval(p))

/*#define R_DT_CIv(p)	R_D_Cval(R_D_qIv(p))		 *  1 - p in qF */
#define R_DT_CIv(p)	(log_p ? (lower_tail ? -expm1(p) : exp(p)) \
			       : R_D_Cval(p))

#define R_DT_exp(x)	R_D_exp(R_D_Lval(x))		/* exp(x) */
#define R_DT_Cexp(x)	R_D_exp(R_D_Cval(x))		/* exp(1 - x) */

#define R_DT_log(p)	(lower_tail? R_D_log(p) : R_D_LExp(p))/* log(p) in qF */
#define R_DT_Clog(p)	(lower_tail? R_D_LExp(p): R_D_log(p))/* log(1-p) in qF*/
#define R_DT_Log(p)	(lower_tail? (p) : R_Log1_Exp(p))
// ==   R_DT_log when we already "know" log_p == TRUE


#define R_Q_P01_check(p)			\
    if ((log_p	&& p > 0) ||			\
	(!log_p && (p < 0 || p > 1)) )		\
	ML_ERR_return_NAN

/* Do the boundaries exactly for q*() functions :
 * Often  _LEFT_ = ML_NEGINF , and very often _RIGHT_ = ML_POSINF;
 *
 * R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)  :<==>
 *
 *     R_Q_P01_check(p);
 *     if (p == R_DT_0) return _LEFT_ ;
 *     if (p == R_DT_1) return _RIGHT_;
 *
 * the following implementation should be more efficient (less tests):
 */
#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)		\
    if (log_p) {					\
	if(p > 0)					\
	    ML_ERR_return_NAN;				\
	if(p == 0) /* upper bound*/			\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
	if(p == ML_NEGINF)				\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
    }							\
    else { /* !log_p */					\
	if(p < 0 || p > 1)				\
	    ML_ERR_return_NAN;				\
	if(p == 0)					\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
	if(p == 1)					\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
    }

#define R_P_bounds_01(x, x_min, x_max) 	\
    if(x <= x_min) return R_DT_0;		\
    if(x >= x_max) return R_DT_1
 /* is typically not quite optimal for (-Inf,Inf) where
  * you'd rather have */
#define R_P_bounds_Inf_01(x)			\
    if(!R_FINITE(x)) {				\
	if (x > 0) return R_DT_1;		\
	/* x < 0 */return R_DT_0;		\
    }



  /* additions for density functions (C.Loader) */
#define R_D_fexp(f,x)     (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))

/* [neg]ative or [non int]eger : */
#define R_D_negInonint(x) (x < 0. || R_nonint(x))

// for discrete d<distr>(x, ...) :
#define R_D_nonint_check(x) 				\
   if(R_nonint(x)) {					\
	MATHLIB_WARNING("non-integer x = %f", x);	\
	return R_D__0;					\
   }
//////////////////////////////////////////////////////////////

double qnorm5(double p)//, double mu, double sigma, int lower_tail, int log_p)
{
    double mu = 0.0;
    double sigma = 1.0;

    int lower_tail = 1.0;
    int log_p = 0;

    double p_, q, r, val;

//#ifdef IEEE_754
//    if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma))
//        return p + mu + sigma;
//#endif
//    R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);
//
//    if (sigma < 0)	ML_ERR_return_NAN;
//    if (sigma == 0)	return mu;
   p_ = R_DT_qIv(p);/* real lower_tail prob. p */
//    log_p ? (lower_tail ? exp(p) : -expm1(p)) : (lower_tail ? (p) : (0.5 - (p)+0.5)); //R_D_Lval(p))


    log_p = p;

    q = p_ - 0.5;

//#ifdef DEBUG_qnorm
//    REprintf("qnorm(p=%10.7g, m=%g, s=%g, l.t.= %d, log= %d): q = %g\n",
//        p, mu, sigma, lower_tail, log_p, q);
//#endif


    ///-- use AS 241 --- 
    /// double ppnd16_(double *p, long *ifault)/
    ///      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

    //        Produces the normal deviate Z corresponding to a given lower
    //        tail area of P; Z is accurate to about 1 part in 10**16.

    //        (original fortran code used PARAMETER(..) for the coefficients
    //         and provided hash codes for checking them...)
    //
    // 0.075 <= p <= 0.925
    if (fabs(q) <= .425) {
        r = .180625 - q * q;
        val =
            q * (((((((r * 2509.0809287301226727 +
                33430.575583588128105) * r + 67265.770927008700853) * r +
                45921.953931549871457) * r + 13731.693765509461125) * r +
                1971.5909503065514427) * r + 133.14166789178437745) * r +
                3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                28729.085735721942674) * r + 39307.89580009271061) * r +
                21213.794301586595867) * r + 5394.1960214247511077) * r +
                687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else { /* closer than 0.075 from {0,1} boundary */

    /* r = min(p, 1-p) < 0.075 */
        if (q > 0)
        {
            //r = R_DT_CIv(p);/* 1-p */

         //   r = (log_p ? (lower_tail ? -expm1(p) : exp(p)) : (lower_tail ? (0.5 - (p)+0.5) : (p));// R_D_Cval(p))
         //   r = (log_p ? (lower_tail ? -expm1(p) : exp(p)) : (lower_tail ? (0.5 - (p)+0.5) : (p));// R_D_Cval(p))
            r = (log_p ? exp(p) : p);

        }
        else
            r = p_;/* = R_DT_Iv(p) ^=  p */

        r = sqrt(-((log_p &&
            ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ?
            p : /* else */ log(r)));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */
//#ifdef DEBUG_qnorm
//        REprintf("\t close to 0 or 1: r = %7g\n", r);
//#endif
        // <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 
        if (r <= 5.) { 
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                .0227238449892691845833) * r + .24178072517745061177) *
                r + 1.27045825245236838258) * r +
                3.64784832476320460504) * r + 5.7694972214606914055) *
                r + 4.6303378461565452959) * r +
                1.42343711074968357734)
                / (((((((r *
                    1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                    r + .0151986665636164571966) * r +
                    .14810397642748007459) * r + .68976733498510000455) *
                    r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.);
        }
        else { /* very close to  0 or 1 */
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                2.71155556874348757815e-5) * r +
                .0012426609473880784386) * r + .026532189526576123093) *
                r + .29656057182850489123) * r +
                1.7848265399172913358) * r + 5.4637849111641143699) *
                r + 6.6579046435011037772)
                / (((((((r *
                    2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
                    r + 1.8463183175100546818e-5) * r +
                    7.868691311456132591e-4) * r + .0148753612908506148525)
                    * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
        }

        if (q < 0.0)
            val = -val;
        /* return (q >= 0.)? r : -r ;*/
    }
    //return mu + sigma * val;
    return val;
}




double qnorm6(double p)
{

   double r, val;   
   double   q = p - 0.5;


    if (fabs(q) <= .425)
    {
        r = .180625 - q * q;
        val =
            q * (((((((r * 2509.0809287301226727 +
                33430.575583588128105) * r + 67265.770927008700853) * r +
                45921.953931549871457) * r + 13731.693765509461125) * r +
                1971.5909503065514427) * r + 133.14166789178437745) * r +
                3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                28729.085735721942674) * r + 39307.89580009271061) * r +
                21213.794301586595867) * r + 5394.1960214247511077) * r +
                687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else 
    { /* closer than 0.075 from {0,1} boundary */

       // r = min(p, 1-p) < 0.075 
        r = std::min(p, 1 - p);
        r = sqrt(-log(r));

        double valMult = (q < 0.0) ? -1.0 : 1.0;

        // <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 
        if (r <= 5.) 
        {
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                .0227238449892691845833) * r + .24178072517745061177) *
                r + 1.27045825245236838258) * r +
                3.64784832476320460504) * r + 5.7694972214606914055) *
                r + 4.6303378461565452959) * r +
                1.42343711074968357734)
                / (((((((r *
                    1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                    r + .0151986665636164571966) * r +
                    .14810397642748007459) * r + .68976733498510000455) *
                    r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.);
        }
        else
        { /* very close to  0 or 1 */
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                2.71155556874348757815e-5) * r +
                .0012426609473880784386) * r + .026532189526576123093) *
                r + .29656057182850489123) * r +
                1.7848265399172913358) * r + 5.4637849111641143699) *
                r + 6.6579046435011037772)
                / (((((((r *
                    2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
                    r + 1.8463183175100546818e-5) * r +
                    7.868691311456132591e-4) * r + .0148753612908506148525)
                    * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
        }

        val *= valMult;
    }

    return val;
}



double qnorm7(double p)
{

    double r, val;
    double   q = p - 0.5;


    static long double a[] = { 2509.0809287301226727 , 33430.575583588128105, 67265.770927008700853, 45921.953931549871457, 13731.693765509461125,  1971.5909503065514427, 133.14166789178437745,3.387132872796366608 };
    static long double b[] = { 5226.495278852854561, 28729.085735721942674,    39307.89580009271061, 21213.794301586595867, 5394.1960214247511077,   687.1870074920579083, 42.313330701600911252 };

    static long double c[] = { 7.7454501427834140764e-4 , .0227238449892691845833 ,.24178072517745061177, 1.27045825245236838258 ,  3.64784832476320460504, 5.7694972214606914055, 4.6303378461565452959, 1.42343711074968357734 };
    static long double d[] = { 1.05075007164441684324e-9 , 5.475938084995344946e-4, .0151986665636164571966, .14810397642748007459, .68976733498510000455,  1.6763848301838038494,  2.05319162663775882187,1. };


    static long double e[] = { 2.01033439929228813265e-7 ,   2.71155556874348757815e-5,   .0012426609473880784386, .026532189526576123093, .29656057182850489123,   1.7848265399172913358, 5.4637849111641143699, 6.6579046435011037772 };

    static long double f[] = { 2.04426310338993978564e-15 , 1.4215117583164458887e-7, 1.8463183175100546818e-5,  7.868691311456132591e-4, .0148753612908506148525,.13692988092273580531, .59983220655588793769, 1. };
        
        



   if (fabs(q) <= .425)
    {
        r = .180625 - q * q;
         val =

         q* (((((((r * a[0] +
             a[1]) * r + a[2]) * r +
             a[3]) * r + a[4]) * r +
             a[5]) * r + a[6]) * r +
             a[7])
             / (((((((r * b[0] +
                 b[1]) * r + b[2]) * r +
                 b[3]) * r +b[4]) * r +
                 b[5]) * r + b[6]) * r + 1.);

    }
    else
    { /* closer than 0.075 from {0,1} boundary */

       // r = min(p, 1-p) < 0.075 
        r = std::min(p, 1 - p);
        r = sqrt(-log(r));

        double valMult = (q < 0.0) ? -1.0 : 1.0;

        // <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 
        if (r <= 5.)
        {
            r += -1.6;
    /*        val = (((((((r * 7.7454501427834140764e-4 +
                .0227238449892691845833) * r + .24178072517745061177) *
                r + 1.27045825245236838258) * r +
                3.64784832476320460504) * r + 5.7694972214606914055) *
                r + 4.6303378461565452959) * r +
                1.42343711074968357734)
                / (((((((r *
                    1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                    r + .0151986665636164571966) * r +
                    .14810397642748007459) * r + .68976733498510000455) *
                    r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.);
    */        

            val = (((((((r * c[0] +  c[1]) * r + c[2]) * r + c[3]) * r +  c[4]) * r + c[5]) *  r + c[6]) * r +    c[7])
                / (((((((r * d[0] + d[1]) *  r + d[2]) * r + d[3]) * r + d[4]) *  r + d[5]) * r +  d[6]) * r + 1.);


        }
        else
        { /* very close to  0 or 1 */
            r += -5.;

            /*
            val = (((((((r * 2.01033439929228813265e-7 +
                2.71155556874348757815e-5) * r +
                .0012426609473880784386) * r + .026532189526576123093) *
                r + .29656057182850489123) * r +
                1.7848265399172913358) * r + 5.4637849111641143699) *
                r + 6.6579046435011037772)
                / (((((((r *
                    2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
                    r + 1.8463183175100546818e-5) * r +
                    7.868691311456132591e-4) * r + .0148753612908506148525)
                    * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
            */

            val = (((((((r * e[0] + e[1]) * r +  e[2]) * r + e[3]) * r + e[4]) * r +  e[5]) * r + e[6]) * r + e[7])
                / (((((((r * f[0] + f[1]) * r + f[2])  * r + f[3]) * r + f[4]) * r +  f[5]) * r + f[6]) * r + 1.);

        }

        val *= valMult;
    }

    return val;
}



/*

double qnorm8(double p)
{
  
    static long double a[] = { 2509.0809287301226727 , 33430.575583588128105, 67265.770927008700853, 45921.953931549871457, 13731.693765509461125,  1971.5909503065514427, 133.14166789178437745,3.387132872796366608 };
    static long double b[] = { 5226.495278852854561, 28729.085735721942674,    39307.89580009271061, 21213.794301586595867, 5394.1960214247511077,   687.1870074920579083, 42.313330701600911252 };
    static long double c[] = { 7.7454501427834140764e-4 , .0227238449892691845833 ,.24178072517745061177, 1.27045825245236838258 ,  3.64784832476320460504, 5.7694972214606914055, 4.6303378461565452959, 1.42343711074968357734 };
    static long double d[] = { 1.05075007164441684324e-9 , 5.475938084995344946e-4, .0151986665636164571966, .14810397642748007459, .68976733498510000455,  1.6763848301838038494,  2.05319162663775882187,1. };
    static long double e[] = { 2.01033439929228813265e-7 ,   2.71155556874348757815e-5,   .0012426609473880784386, .026532189526576123093, .29656057182850489123,   1.7848265399172913358, 5.4637849111641143699, 6.6579046435011037772 };
    static long double f[] = { 2.04426310338993978564e-15 , 1.4215117583164458887e-7, 1.8463183175100546818e-5,  7.868691311456132591e-4, .0148753612908506148525,.13692988092273580531, .59983220655588793769, 1. };


    double val = 0.0;
    double   q = p - 0.5;


    if (fabs(q) <= .425)
    {
        double r = .180625 - q * q;
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
       double  r = std::min(p, 1 - p);
        r = sqrt(-log(r));

        double valMult = (q < 0.0) ? -1.0 : 1.0;

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

        val *= valMult;
    }

    return val;
}

*/