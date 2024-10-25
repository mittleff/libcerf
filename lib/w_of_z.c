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

//--- The following code is generated by devtool/pro_wofz_taylor.py; do not edit
// clang-format off
alignas(64) static const double TaylorCoeffs[64] = {
 0x0.929384f76a0de8p-3, 0x0.8e143a44f9a258p-3, -0x0.8fe9564e0d7190p-5, 0x0.d6d71f181a5048p-10, 0x0.819149ec86c628p-8, -0x0.9c339df5a068d8p-8, 0x0.af4193513c0d28p-13, 0x0.8cdec0fbaa91e0p-10, -0x0.a9abb5b03b5ec0p-13, -0x0.d26aae063555a8p-14, 0x0.84263f053e0400p-15, -0x0.afa55399829e20p-17, -0x0.885284f2f110f0p-19, 0x0.b50489828dee70p-18, -0x0.8a3d375fa09ca8p-21, -0x0.e4e10b4c2563a8p-21, 0x0.b6bac419b127a0p-23, 0x0.a2a5f4d5446ee0p-26, -0x0.abccd634aa9640p-26, 0x0.b79920a9369d98p-26, -0x0.80c74446749290p-30, -0x0.a7df4bb92a09d8p-28, 0x0.d22aced7f291b0p-31, 0x0.bfa2bbb51ddbb0p-32, -0x0.853a279b18f178p-33, 0x0.9eb78ea9b86b98p-34, 0x0.be80076341b838p-40, -0x0.ce84e0640acfa8p-36, 0x0.db8faf5b322390p-39, 0x0.a4184eee414ba0p-39, -0x0.a920cfad8aca98p-41, 0x0.fab9aee345dee0p-43, 0x0.c3f8c0b2a1ec10p-46, -0x0.e697abe61fe2f0p-44, 0x0.dce0fd4c6dfd28p-47, 0x0.d3b71f2c894d28p-47, -0x0.be7ce0d77bd1d8p-49, 0x0.cb4fdec2036eb0p-51, 0x0.89ae30b9a65330p-53, -0x0.f2b69b053dc6c0p-52, 0x0.de26942948a458p-55, 0x0.e7183db14fac28p-55, -0x0.c42b78ad237f00p-57, 0x0.c426f4f31fed60p-59, 0x0.862695d527c5d8p-61, -0x0.f4417aeb8075d8p-60, 0x0.e16a364e76aa60p-63, 0x0.dc1dcc9a0cd418p-63, -0x0.b9eeaa0cc22580p-65, 0x0.e53a64f7f24670p-67, 0x0.b685a3b0f64520p-70, -0x0.ea37a350128fd0p-68, 0x0.e280b0e56cb608p-71, 0x0.b6163ad8449200p-71, -0x0.a0cd5eb1593738p-73, 0x0.8ddaea64a0cfa8p-74, 0x0.bf299208b39340p-80, -0x0.d2c54c39ad2a90p-76, 0x0.da3b3dddbde2a0p-79, 0x0.fb71ce25896890p-80, -0x0.f8545791326800p-82, 0x0.a5853ab0bb2f68p-82, -0x0.babd35783306c0p-87, -0x0.ae134ba2234a50p-84,};
// clang-format on
//--- End of autogenerated code

_cerf_cmplx taylor(double xa, double ya, int jblock)
{
    static const int N = 31;
    /*
    static const int jblock0 = -1;
    static const int L = 3;
    static const int nT = 2*(N+1) + 6; // number of entries per block
    static const int nB = 3 * (1 << (2*(L-1))); // number of blocks per call (i.e. per value of jblock)

    int jx = (1 << L) - 1 - (int)(xa * (1<<(L-jblock)));
    int jy = (1 << L) - 1 - (int)(ya * (1<<(L-jblock)));

    int kblock = ((jx&1) + ((jx&2)<<1) + ((jx&4)<<2) + ((jx&8)<<3)
		  + ((jy&1)<<1) + ((jy&2)<<2) + ((jy&4)<<3) + ((jy&8)<<4)
		  + nB * (jblock - jblock0)) * nT;

    assert(kblock <= (5*nB-1)*nT);
    */
    int kblock = 0;
    const double* T = &TaylorCoeffs[kblock];
    /*
    if (xa<T[nT-4] || xa>T[nT-3] || ya<T[nT-2] || ya>T[nT-1]) {
	printf("jx=%i jy=%i  cx=%i cy=%i  jb=%i kb=%i\n",
	       jx, jy, (jx&1) + ((jx&2)<<1) + ((jx&4)<<2) + ((jx&8)<<3),
	       ((jy&1)<<1) + ((jy&2)<<2) + ((jy&4)<<3) + ((jy&8)<<4), jblock, kblock);
	printf("x %g in %g %g\n", xa, T[nT-4], T[nT-3]);
	printf("y %g in %g %g\n", ya, T[nT-2], T[nT-1]);
	exit(-1);
    }
    */
    _cerf_cmplx ret = C(T[2*N], T[2*N+1]);
    //_cerf_cmplx dz = C(xa, ya) - C(T[nT-6], T[nT-5]);
    _cerf_cmplx dz = C(xa, ya) - C(4., 4.);
    for (int k = N-1; k >= 0; --k)
	ret = ret * dz + C(T[2*k], T[2*k+1]);
    SET_INFO(900, kblock);
    return ret;
}

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

    const double xs = y < 0 ? -creal(z) : creal(z); // compute for -z if y < 0

// ------------------------------------------------------------------------------
// Taylor in y                                                         [ALGO 7??]
// ------------------------------------------------------------------------------

    if (ya < .25 + xa/8) {
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
// Taylor around some z                                                [ALGO 9??]
// ------------------------------------------------------------------------------

    ret = taylor(xa, ya, 0);

    if (y < 0) {
	SET_ALGO(cerf_algorithm + 1);
	return 2.0 * cexp(C((ya - xs) * (xs + ya), 2*xs*y)) - ret;
    } else
	return ret;

// ------------------------------------------------------------------------------
// Failed                                                                [ALGO 0]
// ------------------------------------------------------------------------------

    SET_INFO(0, 0);
    return 0;

} // w_of_z
