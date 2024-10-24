/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File w_of_z.c:
 *   Computation of Faddeeva's complex scaled error function,
 *      w(z) = exp(-z^2) * erfc(-i*z),
 *   nameless function (7.1.3) of Abramowitz&Stegun (1964),
 *   also known as the plasma dispersion function.
 *
 *   This implementation uses a combination of different algorithms.
 *   See man 3 w_of_z for references.
 *
 * Copyright:
 *   (C) 2012 Massachusetts Institute of Technology
 *   (C) 2013 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   Permission is hereby granted, free of charge, to any person obtaining
 *   a copy of this software and associated documentation files (the
 *   "Software"), to deal in the Software without restriction, including
 *   without limitation the rights to use, copy, modify, merge, publish,
 *   distribute, sublicense, and/or sell copies of the Software, and to
 *   permit persons to whom the Software is furnished to do so, subject to
 *   the following conditions:
 *
 *   The above copyright notice and this permission notice shall be
 *   included in all copies or substantial portions of the Software.
 *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *   LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *   OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * Authors:
 *   Steven G. Johnson, Massachusetts Institute of Technology, 2012
 *   Joachim Wuttke, Forschungszentrum Jülich, 2013, 2024
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/libcerf
 *
 * Revision history:
 *   ../CHANGELOG
 *
 * Man page:
 *   w_of_z(3)
 */

/*
  Computes various error functions (erf, erfc, erfi, erfcx),
  including the Dawson integral, in the complex plane, based
  on algorithms for the computation of the Faddeeva function
  w(z) = exp(-z^2) * erfc(-i*z).
  Given w(z), the error functions are mostly straightforward
  to compute, except for certain regions where we have to
  switch to Taylor expansions to avoid cancellation errors
  [e.g. near the origin for erf(z)].

*/

#include "cerf.h"
#include "defs.h" // defines _cerf_cmplx, NaN, C, cexp, ...
#include <float.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#ifdef CERF_INTROSPECT
EXPORT int cerf_algorithm;
EXPORT int cerf_nofterms;
#endif

/******************************************************************************/
/*  auxiliary functions                                                       */
/******************************************************************************/

static inline double sinc(double x, double sinx) {
    // return sinc(x) = sin(x)/x, given both x and sin(x)
    // [since we only use this in cases where sin(x) has already been computed]
    return fabs(x) < 1e-4 ? 1 - (0.1666666666666666666667) * x * x : sinx / x;
}

/******************************************************************************/
/* precomputed table of expa2n2[n-1] = exp(-a^2 * n^2)                        */
/******************************************************************************/

static const double expa2n2[] = {
    7.78800783071404878e-1,
    3.67879441171442334e-1,
    1.05399224561864333e-1,
    1.83156388887341787e-2,
    1.9304541362277093e-3,
    1.23409804086679561e-4,
    4.78511739212900875e-6,
    1.12535174719259116e-7,
    1.60522805518561165e-9,
    1.38879438649640209e-11,
    7.28772409581969219e-14,
    2.31952283024356963e-16,
    4.4777324417183015e-19,
    5.2428856633634639e-22,
    3.72336312175051061e-25,
    1.60381089054863793e-28,
    4.19009319449439736e-32,
    6.63967719958073481e-36,
    6.38150344806079078e-40,
    3.72007597602083612e-44,
    1.31532589485746445e-48,
    2.8207700884601352e-53,
    3.66905961542916342e-58,
    2.89464031164830029e-63,
    1.38511936992260166e-68,
    4.02006021574335523e-74,
    7.07669817542954913e-80,
    7.55581901971196089e-86,
    4.8931122620973616e-92,
    1.92194772782384913e-98,
    4.57878996929157753e-105,
    6.61626105670948528e-112,
    5.79865541856403774e-119,
    3.08244069694909812e-126,
    9.93836441348368381e-134,
    1.9435148500492928e-141,
    2.30522631523556484e-149,
    1.65841047768114525e-157,
    7.23641151924800958e-166,
    1.91516959671400568e-174,
    3.07428396708383541e-183,
    2.99318445226019286e-192,
    1.76756641155092425e-201,
    6.33097733621059151e-211,
    1.37536679932640647e-220,
    1.81225402579399229e-230,
    1.44834614899898877e-240,
    7.02066779850473475e-251,
    2.06413091092950948e-261,
    3.68085585480180036e-272,
    3.9811921806329143e-283,
    2.61174176128405546e-294,
    1.03920226214308252e-305,
    0.0 // underflow (will force termination of loops, and thereby prevent reading past array end)
}; // expa2n2

/******************************************************************************/
/*  w_of_z, Faddeeva's scaled complex error function                          */
/******************************************************************************/

_cerf_cmplx w_of_z(_cerf_cmplx z) {
    SET_INFO(-1, -1);

    const double x = creal(z);
    const double xa = fabs(x);
    const double y = cimag(z);
    const double ya = fabs(y);
    const double z2 = xa*xa + y*y;

    const double ispi = 0.5641895835477562869; // 1 / sqrt(pi)

// ------------------------------------------------------------------------------
// Case |y| << |x|                                                     [ALGO 3??]
// ------------------------------------------------------------------------------

//   If |y| << |x|, we get |Re w| << |Im w|, which implies that an algorithm
//   can be very accurate in terms of the complex norm, with relative error
//   computed as |w - w_0| / |w_0|, and still Re w can be wrong by orders of
//   magnitude. Therefore this special case is treated beforehand.

//   We start from the real axis where w(x) = exp(-x^2) + i*im_w_of_x(x).
//   To compute w(x+iy) we assume that |y| is negligible when added to |x|.
//   We obtain Re w = exp(-x^2) + 2 y (x Im w(x) - 1 / sqrt(pi)).

    if (ya < 1e-16 * xa) {
	const double wi = im_w_of_x(x);
	SET_ALGO(cerf_algorithm + 300);
        const double e2 = xa > 27. ? 0. : exp(-xa*xa); // prevent underflow
	if (ya == 0)
	    return C(e2, wi); // also works for x=+-inf
	return C(e2 + y*(2*(x*wi - ispi)), wi);
    }

// ------------------------------------------------------------------------------
// Case |x| << |y|                                                     [ALGO 4??]
// ------------------------------------------------------------------------------

//   If |x| << |y|, we get |Im w| << |Re w|, which implies that an algorithm
//   can be very accurate in terms of the complex norm, with relative error
//   computed as |w - w_0| / |w_0|, and still Im w can be wrong by orders of
//   magnitude. Therefore this special case is treated beforehand.

//   We start from the imaginary axis where w(iy) = erfcx(y).
//   To compute w(iy+x) we assume that |x| is negligible when added to |y|.
//   We obtain Im w = 2 x (1 / sqrt(pi) - y erfcx(y)).

    if (xa < 1e-16 * ya) {
	const double wr = erfcx(y);
	SET_ALGO(cerf_algorithm + 400);
	if (xa == 0)
	    return C(wr, 0); // also works for y=+inf
	return C(wr, x*(2*(ispi - y*wr)));
    }

// ------------------------------------------------------------------------------
// Case |z| -> 0: Maclaurin series                                     [ALGO 210]
// ------------------------------------------------------------------------------

    if (z2 < .26) {
        if (z2 < .00689) {
            if (z2 < 4e-7) {
		SET_INFO(210, 5);
                return ((((
                              + C(+5.0000000000000000e-01, 0) ) * z // z^4
                          + C(0, -7.5225277806367508e-01) ) * z // z^3
                         + C(-1.0000000000000000e+00, 0) ) * z // z^2
                        + C(0, +1.1283791670955126e+00) ) * z // z^1
                    + 1.;
            }

	    SET_INFO(210, 14);
            return (((((((((((((
                                   + C(0, +5.3440090793734269e-04) ) * z // z^13
                               + C(+1.3888888888888889e-03, 0) ) * z // z^12
                              + C(0, -3.4736059015927274e-03) ) * z // z^11
                             + C(-8.3333333333333332e-03, 0) ) * z // z^10
                            + C(0, +1.9104832458760001e-02) ) * z // z^9
                           + C(+4.1666666666666664e-02, 0) ) * z // z^8
                          + C(0, -8.5971746064419999e-02) ) * z // z^7
                         + C(-1.6666666666666666e-01, 0) ) * z // z^6
                        + C(0, +3.0090111122547003e-01) ) * z // z^5
                       + C(+5.0000000000000000e-01, 0) ) * z // z^4
                      + C(0, -7.5225277806367508e-01) ) * z // z^3
                     + C(-1.0000000000000000e+00, 0) ) * z // z^2
                    + C(0, +1.1283791670955126e+00) ) * z // z^1
                + 1.;
        }

        if (z2 < .074) {
	    SET_INFO(210, 20);
            return (((((((((((((((((((
                                         + C(0, -8.8239572002038009e-07) ) * z // z^19
                                     + C(-2.7557319223985893e-06, 0) ) * z // z^18
                                    + C(0, +8.3827593401936105e-06) ) * z // z^17
                                   + C(+2.4801587301587302e-05, 0) ) * z // z^16
                                  + C(0, -7.1253454391645692e-05) ) * z // z^15
                                 + C(-1.9841269841269841e-04, 0) ) * z // z^14
                                + C(0, +5.3440090793734269e-04) ) * z // z^13
                               + C(+1.3888888888888889e-03, 0) ) * z // z^12
                              + C(0, -3.4736059015927274e-03) ) * z // z^11
                             + C(-8.3333333333333332e-03, 0) ) * z // z^10
                            + C(0, +1.9104832458760001e-02) ) * z // z^9
                           + C(+4.1666666666666664e-02, 0) ) * z // z^8
                          + C(0, -8.5971746064419999e-02) ) * z // z^7
                         + C(-1.6666666666666666e-01, 0) ) * z // z^6
                        + C(0, +3.0090111122547003e-01) ) * z // z^5
                       + C(+5.0000000000000000e-01, 0) ) * z // z^4
                      + C(0, -7.5225277806367508e-01) ) * z // z^3
                     + C(-1.0000000000000000e+00, 0) ) * z // z^2
                    + C(0, +1.1283791670955126e+00) ) * z // z^1
                + 1.;
        }

	SET_INFO(210, 26);
        return (((((((((((((((((((((((((
                                           + C(0, +5.8461000084165970e-10) ) * z // z^25
                                       + C(+2.0876756987868100e-09, 0) ) * z // z^24
                                      + C(0, -7.3076250105207460e-09) ) * z // z^23
                                     + C(-2.5052108385441720e-08, 0) ) * z // z^22
                                    + C(0, +8.4037687620988577e-08) ) * z // z^21
                                   + C(+2.7557319223985888e-07, 0) ) * z // z^20
                                  + C(0, -8.8239572002038009e-07) ) * z // z^19
                                 + C(-2.7557319223985893e-06, 0) ) * z // z^18
                                + C(0, +8.3827593401936105e-06) ) * z // z^17
                               + C(+2.4801587301587302e-05, 0) ) * z // z^16
                              + C(0, -7.1253454391645692e-05) ) * z // z^15
                             + C(-1.9841269841269841e-04, 0) ) * z // z^14
                            + C(0, +5.3440090793734269e-04) ) * z // z^13
                           + C(+1.3888888888888889e-03, 0) ) * z // z^12
                          + C(0, -3.4736059015927274e-03) ) * z // z^11
                         + C(-8.3333333333333332e-03, 0) ) * z // z^10
                        + C(0, +1.9104832458760001e-02) ) * z // z^9
                       + C(+4.1666666666666664e-02, 0) ) * z // z^8
                      + C(0, -8.5971746064419999e-02) ) * z // z^7
                     + C(-1.6666666666666666e-01, 0) ) * z // z^6
                    + C(0, +3.0090111122547003e-01) ) * z // z^5
                   + C(+5.0000000000000000e-01, 0) ) * z // z^4
                  + C(0, -7.5225277806367508e-01) ) * z // z^3
                 + C(-1.0000000000000000e+00, 0) ) * z // z^2
                + C(0, +1.1283791670955126e+00) ) * z // z^1
            + 1.;
    }

    _cerf_cmplx ret = 0.; // return value

// ------------------------------------------------------------------------------
// Case |z| -> infty: Asymptotic expansion                        [ALGO 100, 22?]
// ------------------------------------------------------------------------------

    if (z2 > 49) {
	const double xs = y < 0 ? -creal(z) : creal(z); // compute for -z if y < 0

	if (z2 > 4.8e15) {
	    // Scale to prevent overflow.
	    if (xa > ya) {
		SET_INFO(222, 1);
		const double yax = ya / xs;
		const double denom = 0.56418958354775629 / (xs + yax*ya);
		ret = C(denom*yax, denom);
	    } else if (isinf(ya)) {
		SET_INFO(100, 1);
		return ((isnan(xa) || y < 0) ? C(NaN, NaN) : C(0, 0));
	    } else {
		SET_INFO(224, 1);
		const double xya = xs / ya;
		const double denom = 0.56418958354775629 / (xya*xs + ya);
		ret = C(denom, denom*xya);
	    }

	} else {
	    const double zm2 = 1 / z2;                                   // 1/|z|^2
	    const _cerf_cmplx r = C(ya*zm2, xs*zm2);                     // i/z
	    const double zm4 = zm2 * zm2;                                // 1/|z|^4
	    const _cerf_cmplx r2 = C(zm4*(xs+ya)*(xs-ya), -2*zm4*xs*ya); // 1/z^2

            if (z2 > 540) {
		if (z2 > 22500) {
		    SET_INFO(220, 4);
		    ret = ((((
				 + 1.0578554691520430e+00) * r2 // n=3
			     + 4.2314218766081724e-01) * r2 // n=2
			    + 2.8209479177387814e-01) * r2 // n=1
			   + 5.6418958354775628e-01) * r; // n=0

		} else {
		    SET_INFO(220, 12);
		    ret = ((((((((((((
					 + 3.7877040075087948e+06) * r2 // n=11
				     + 3.6073371500083758e+05) * r2 // n=10
				    + 3.7971970000088164e+04) * r2 // n=9
				   + 4.4672905882456671e+03) * r2 // n=8
				  + 5.9563874509942218e+02) * r2 // n=7
				 + 9.1636730015295726e+01) * r2 // n=6
				+ 1.6661223639144676e+01) * r2 // n=5
			       + 3.7024941420321507e+00) * r2 // n=4
			      + 1.0578554691520430e+00) * r2 // n=3
			     + 4.2314218766081724e-01) * r2 // n=2
			    + 2.8209479177387814e-01) * r2 // n=1
			   + 5.6418958354775628e-01) * r; // n=0
		}

            } else {
		SET_INFO(220, 20);
                ret = ((((((((((((((((((((
                                             + 8.8249260943025370e+15) * r2 // n=19
                                         + 4.7702303212446150e+14) * r2 // n=18
                                        + 2.7258458978540656e+13) * r2 // n=17
                                       + 1.6520278168812520e+12) * r2 // n=16
                                      + 1.0658243979879044e+11) * r2 // n=15
                                     + 7.3505130895717545e+09) * r2 // n=14
                                    + 5.4448245107938921e+08) * r2 // n=13
                                   + 4.3558596086351141e+07) * r2 // n=12
                                  + 3.7877040075087948e+06) * r2 // n=11
                                 + 3.6073371500083758e+05) * r2 // n=10
                                + 3.7971970000088164e+04) * r2 // n=9
                               + 4.4672905882456671e+03) * r2 // n=8
                              + 5.9563874509942218e+02) * r2 // n=7
                             + 9.1636730015295726e+01) * r2 // n=6
                            + 1.6661223639144676e+01) * r2 // n=5
                           + 3.7024941420321507e+00) * r2 // n=4
                          + 1.0578554691520430e+00) * r2 // n=3
                         + 4.2314218766081724e-01) * r2 // n=2
                        + 2.8209479177387814e-01) * r2 // n=1
                       + 5.6418958354775628e-01) * r; // n=0
            }
        }
        if (y < 0) {
            // Use w(z) = 2.0*exp(-z*z) - w(-z),
            // but be careful of overflow in exp(-z*z) = exp(-(xs*xs-ya*ya) -2*i*xs*ya)
	    SET_ALGO(cerf_algorithm + 1);
            return 2.0 * cexp(C((ya - xs) * (xs + ya), 2*xs*y)) - ret;
        } else
            return ret;
    }

    const double relerr = DBL_EPSILON;
    const double a = 0.5;  // smaller than pi / sqrt(-log(eps*0.5))
    const double a2_pi = 0.31830988618379067;  // (2/pi) * a;

    double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0;

// ------------------------------------------------------------------------------
// Taylor in y                                                         [ALGO 7??]
// ------------------------------------------------------------------------------

    if (ya < .25 + xa/8) {
        const double xs = y < 0 ? -creal(z) : creal(z); // compute for -z if y < 0
	const double wi = im_w_of_x(xs);
	SET_INFO(cerf_algorithm + 700, 0);
        const double e2 = xa > 27. ? 0. : exp(-xa*xa); // prevent underflow
	ret = C(e2 + y*(2*(x*wi - ispi)), wi);
	static const int NW = 40;
	_cerf_cmplx W[NW];
	W[0] = C(e2, wi);
	W[1] = -2 * xs * W[0] + C(0, 2*ispi);
	for (int n = 2; n < NW; ++n)
	    W[n] = -2. * (W[n-1]*xs + W[n-2]*(n-1));
	ret = 0;
	_cerf_cmplx h = C(0, ya);
	for (int n = NW - 1; n >= 0; --n)
	    ret = ret * h / (n+1.) + W[n];

        if (y < 0) {
            // Use w(z) = 2.0*exp(-z*z) - w(-z),
            // but be careful of overflow in exp(-z*z) = exp(-(xs*xs-ya*ya) -2*i*xs*ya)
	    SET_ALGO(cerf_algorithm + 1);
            return 2.0 * cexp(C((ya - xs) * (xs + ya), 2*xs*y)) - ret;
        } else
            return ret;
    }

// ------------------------------------------------------------------------------
// Taylor in x                                                         [ALGO 8??]
// ------------------------------------------------------------------------------

    if (xa < .25 + ya/8) {
        const double xs = y < 0 ? -creal(z) : creal(z); // compute for -z if y < 0
	const double wr = erfcx(ya);
	SET_INFO(cerf_algorithm + 800, 0);
	static const int NW = 40;
	_cerf_cmplx z0 = C(0., ya);
	_cerf_cmplx W[NW];
	W[0] = wr;
	W[1] = -2. * z0 * W[0] + C(0, 2*ispi);
	for (int n = 2; n < NW; ++n)
	    W[n] = -2. * (W[n-1]*z0 + W[n-2]*(n-1));
	ret = 0;
	for (int n = NW - 1; n >= 0; --n)
	    ret = ret * xs / (n+1.) + W[n];

        if (y < 0) {
            // Use w(z) = 2.0*exp(-z*z) - w(-z),
            // but be careful of overflow in exp(-z*z) = exp(-(xs*xs-ya*ya) -2*i*xs*ya)
	    SET_ALGO(cerf_algorithm + 1);
            return 2.0 * cexp(C((ya - xs) * (xs + ya), 2*xs*y)) - ret;
        } else
            return ret;
    }

// ------------------------------------------------------------------------------
// Intermediate case: continued fraction expansion                     [ALGO 230]
// ------------------------------------------------------------------------------

// Continued-fraction expansion,
// similar to those described by Gautschi (1970) and Poppe & Wijers (1990).
//
// Prefered for large z because it is fast.
//
// However, as pointed out by M. Zaghloul, the continued fraction seems to
// give a large relative error in Re w(z) for |x| ~ 6 and small |y|.
// In this region, excluded by the additional conditions in the if clause above,
// we fall back to ACM algorithm 916 below.
//
// Poppe & Wijers suggest using a number of terms
// nu = 3 + 1442 / (26*rho + 77)
// where rho = sqrt((x/x0)^2 + (y/y0)^2) where x0=6.3, y0=4.4.
// (They only use this expansion for rho >= 1, but rho a little less
// than 1 seems okay too.)
// Instead, I [SGJ] did my own fit to a slightly different function
// that avoids the hypotenuse calculation, using NLopt to minimize
// the sum of the squares of the errors in nu with the constraint
// that the estimated nu be >= minimum nu to attain machine precision.
// I also separate the regions where nu == 2 and nu == 1.

    if (xa > 6 && ya > 0.1) {

        const double xs = y < 0 ? -creal(z) : creal(z); // compute for -z if y < 0

        if (xa + ya > 4000) {                      // nu <= 2
	    assert(0); // must never be reached, has been replaced by the above

        } else {
//	    // Forward recursion for numerator and denominator (GiST07 ch 6.2)
//	    SET_INFO(230, 0);
//	    _cerf_cmplx A2 = C(0., 0.);
//	    _cerf_cmplx B2 = C(1., 0.);
//	    _cerf_cmplx A1 = C(-.5, 0.);
//	    _cerf_cmplx B1 = z;
//	    _cerf_cmplx s = 1/z/2;
//	    for (int n = 2; n < 130; ++n) {
//		_cerf_cmplx A = A1 - A2*n*s;
//		_cerf_cmplx B = B1 - B2*n*s;
//		A2 = A1;
//		B2 = B1;
//		A1 = A;
//		B1 = B;
//		printf("%2i %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n",
//		       n, creal(A1), cimag(A1), creal(B1), cimag(B1), creal(A1/B1), cimag(A1/B1));
//	    }
//	    ret = C(0, ispi) / (z + A1/B1);

	    // Modified Lentz algorithm according to Gil Segura Temma 2007, Algorithm 6.2
	    SET_INFO(230, 0);
	    _cerf_cmplx Cn = C(0, 1e-30);
	    _cerf_cmplx Dn = C(0, 0);
	    _cerf_cmplx En = Cn;
	    for (int n = 1; n < 100; ++n) {
		Dn = z - n / 2. * Dn;
		if (cabs(Dn) < 1e-30)
		    Dn = C(0, 0);
		En = z - n / 2. * En;
		if (cabs(En) < 1e-30)
		    En = C(0, 0);
		Dn = 1 / Dn;
		_cerf_cmplx Hn = En * Dn;
		Cn *= Hn;
		if (cabs(Hn-1) < 1e-16) {
		    SET_NTER(n);
		    break;
		}
	    }
	    ret = C(0, ispi) / (z + Cn);
        }

        if (y < 0) {
            // Use w(z) = 2.0*exp(-z*z) - w(-z),
            // but be careful of overflow in exp(-z*z) = exp(-(xs*xs-ya*ya) -2*i*xs*ya)
	    SET_ALGO(cerf_algorithm + 1);
            return 2.0 * cexp(C((ya - xs) * (xs + ya), 2*xs*y)) - ret;
        } else
            return ret;
    }

// ------------------------------------------------------------------------------
// NaN                                                             [ALGO 105 106]
// ------------------------------------------------------------------------------

    if (isnan(xa)) {
	SET_INFO(105, 1);
	return C(xa, xa);
    }
    if (isnan(y)) {
	SET_INFO(106, 1);
	return C(y, y);
    }

// ------------------------------------------------------------------------------
// Intermediate case: ACM algorithm 916 by Zaghloul & Ali          [ALGO 5?? 6??]
// ------------------------------------------------------------------------------

//  ACM algorithm 916 by Zaghloul & Ali (2011), which is generally competitive
//  at small |z|, and more accurate than the Poppe & Wijers expansion in some
//  regions, e.g. in the vicinity of z=1+i.
//
//  Applicability range: The paper seems to suggest x < sqrt(-log(DBL_MIN)),
//  about 26.6, since otherwise exp(-x^2) underflows to zero and sum1,sum2,sum4
//  are zero.  However, long before this occurs, the sum1,sum2,sum4 contributions
//  are negligible in double precision; I [SGJ] find that this happens for
//  x > about 6, for all y.  On the other hand, I find that the case
//  where we compute all of the sums is faster (at least with the
//  precomputed expa2n2 table) until about x=10.

    double prod2ax = 1, prodm2ax = 1;
    double expx2;

    SET_ALGO(0);
    double e2y;

    // x > 5e-4, compute sum4 and sum5 separately
    expx2 = exp(-xa*xa);
    e2y = y < -6  ? 2*exp(y*y-xa*xa) : expx2*erfcx(y); // may SET_ALGO
    const double exp2ax = exp((2*a) * xa);
    const double expm2ax = 1 / exp2ax;
    for (int n = 1;; ++n) {
	const double coef = expa2n2[n - 1] * expx2 / ((a*a) * (n*n) + y*y);
	prod2ax *= exp2ax;
	prodm2ax *= expm2ax;
	sum1 += coef;
	sum2 += coef * prodm2ax;
	sum3 += coef * prod2ax;
	sum4 += (coef * prodm2ax) * (a*n);
	sum5 += (coef * prod2ax) * (a*n);
	// test convergence via sum5, since this sum has the slowest decay;
	// for termination rely on coef[n_max] = 0
	if ((coef * prod2ax) * (a*n) < relerr * sum5) {
	    SET_INFO(500+cerf_algorithm, n);
	    break;
	}
    }

// The second case has the exact expression.
// In the first case, to avoid spurious overflow for large negative y,
// we approximate erfcx(y) by 2*exp(y^2), which is accurate to double precision.
//
// TODO: check exact location of cross-over

    if (y > 5) { // imaginary terms cancel
	SET_ALGO(cerf_algorithm + 1);
	const double sinxy = sin(xa*y);
	ret = (e2y - a2_pi*y*sum1) * cos(2*xa*y) + (a2_pi*xa*expx2) * sinxy * sinc(xa*y, sinxy);
    } else {
	double xs = creal(z);
	const double sinxy = sin(xs*y);
	const double sin2xy = sin(2*xs*y);
	const double cos2xy = cos(2*xs*y);
	const double coef1 = e2y - a2_pi*y*sum1;
	const double coef2 = a2_pi*xs*expx2;
	ret = C(coef1*cos2xy + coef2*sinxy*sinc(xs*y, sinxy),
		coef2*sinc(2*xs*y, sin2xy) - coef1*sin2xy);
    }
    return ret + C((a2_pi/2) * y * (sum2 + sum3), (a2_pi/2) * copysign(sum5 - sum4, creal(z)));

} // w_of_z
