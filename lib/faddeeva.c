/* Library libcerf:
 *   compute complex error functions,
 *   along with Dawson, Faddeeva and Voigt functions
 *
 * File faddeeva.c:
 *   Computation of Faddeeva, Dawson, Voigt, and error functions.
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
 *   Steven G. Johnson, Massachusetts Institute of Technology, 2012, core author
 *   Joachim Wuttke, Forschungszentrum Jülich, 2013, package maintainer
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/libcerf
 *
 * Revision history:
 *   ../CHANGELOG
 *
 * More information:
 *   man 3 faddeeva
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

    -- Note that Algorithm 916 assumes that the erfc(x) function, 
       or rather the scaled function erfcx(x) = exp(x*x)*erfc(x),
       is supplied for REAL arguments x.   I originally used an
       erfcx routine derived from DERFC in SLATEC, but I have
       since replaced it with a much faster routine written by
       me which uses a combination of continued-fraction expansions
       and a lookup table of Chebyshev polynomials.  For speed,
       I implemented a similar algorithm for Im[w(x)] of real x,
       since this comes up frequently in the other error functions.

   If HAVE_CONFIG_H is #defined (e.g. by compiling with -DHAVE_CONFIG_H),
   then we #include "config.h", which is assumed to be a GNU autoconf-style
   header defining HAVE_* macros to indicate the presence of features. In
   particular, if HAVE_ISNAN and HAVE_ISINF are #defined, we use those
   functions in math.h instead of defining our own, and if HAVE_ERF and/or
   HAVE_ERFC are defined we use those functions from <cmath> for erf and
   erfc of real arguments, respectively, instead of defining our own.
*/

#include "cerf.h"

#define _GNU_SOURCE // enable GNU libc NAN extension if possible

#include <float.h>
#include <math.h>

#include "defs.h" // defines cmplx, CMPLX, NaN

/******************************************************************************/
static inline cmplx cpolar(double r, double t)
{
    if (r == 0.0 && !isnan(t))
        return 0.0;
    else
        return C(r * cos(t), r * sin(t));
}

/////////////////////////////////////////////////////////////////////////
// Auxiliary routines to compute other special functions based on w(z)

/******************************************************************************/
// compute erfcx(z) = exp(z^2) erfz(z)
cmplx faddeeva_erfcx(cmplx z, double relerr)
{
    return faddeeva(C(-cimag(z), creal(z)), relerr);
}

/******************************************************************************/
// compute the error function erf(x)
double faddeeva_erf_re(double x)
{
    return erf(x); // C99 supplies erf in math.h
}

/******************************************************************************/
// compute the error function erf(z)
cmplx faddeeva_erf(cmplx z, double relerr)
{
    double x = creal(z), y = cimag(z);

    if (y == 0)
        return C(faddeeva_erf_re(x),
                 y); // preserve sign of 0
    if (x == 0) // handle separately for speed & handling of y = Inf or NaN
        return C(x, // preserve sign of 0
                 /* handle y -> Inf limit manually, since
                    exp(y^2) -> Inf but Im[w(y)] -> 0, so
                    IEEE will give us a NaN when it should be Inf */
                 y*y > 720 ? (y > 0 ? Inf : -Inf)
                 : exp(y*y) * im_w_of_z(y));
  
    double mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
    double mIm_z2 = -2*x*y; // Im(-z^2)
    if (mRe_z2 < -750) // underflow
        return (x >= 0 ? 1.0 : -1.0);

    /* Handle positive and negative x via different formulas,
       using the mirror symmetries of w, to avoid overflow/underflow
       problems from multiplying exponentially large and small quantities. */
    if (x >= 0) {
        if (x < 8e-2) {
            if (fabs(y) < 1e-2)
                goto taylor;
            else if (fabs(mIm_z2) < 5e-3 && x < 5e-3)
                goto taylor_erfi;
        }
        /* don't use complex exp function, since that will produce spurious NaN
           values when multiplying w in an overflow situation. */
        return 1.0 - exp(mRe_z2) *
            (C(cos(mIm_z2), sin(mIm_z2))
             * faddeeva(C(-y,x), relerr));
    }
    else { // x < 0
        if (x > -8e-2) { // duplicate from above to avoid fabs(x) call
            if (fabs(y) < 1e-2)
                goto taylor;
            else if (fabs(mIm_z2) < 5e-3 && x > -5e-3)
                goto taylor_erfi;
        }
        else if (isnan(x))
            return C(NaN, y == 0 ? 0 : NaN);
        /* don't use complex exp function, since that will produce spurious NaN
           values when multiplying w in an overflow situation. */
        return exp(mRe_z2) *
            (C(cos(mIm_z2), sin(mIm_z2))
             * faddeeva(C(y,-x), relerr)) - 1.0;
    }

    // Use Taylor series for small |z|, to avoid cancellation inaccuracy
    //   erf(z) = 2/sqrt(pi) * z * (1 - z^2/3 + z^4/10 - z^6/42 + z^8/216 + ...)
taylor:
    {
        cmplx mz2 = C(mRe_z2, mIm_z2); // -z^2
        return z * (1.1283791670955125739
                    + mz2 * (0.37612638903183752464
                             + mz2 * (0.11283791670955125739
                                      + mz2 * (0.026866170645131251760
                                               + mz2 * 0.0052239776254421878422))));
    }

    /* for small |x| and small |xy|, 
       use Taylor series to avoid cancellation inaccuracy:
       erf(x+iy) = erf(iy)
       + 2*exp(y^2)/sqrt(pi) *
       [ x * (1 - x^2 * (1+2y^2)/3 + x^4 * (3+12y^2+4y^4)/30 + ... 
       - i * x^2 * y * (1 - x^2 * (3+2y^2)/6 + ...) ]
       where:
       erf(iy) = exp(y^2) * Im[w(y)]
    */
taylor_erfi:
    {
        double x2 = x*x, y2 = y*y;
        double expy2 = exp(y2);
        return C
            (expy2 * x * (1.1283791670955125739
                          - x2 * (0.37612638903183752464
                                  + 0.75225277806367504925*y2)
                          + x2*x2 * (0.11283791670955125739
                                     + y2 * (0.45135166683820502956
                                             + 0.15045055561273500986*y2))),
             expy2 * (im_w_of_z(y)
                      - x2*y * (1.1283791670955125739 
                                - x2 * (0.56418958354775628695 
                                        + 0.37612638903183752464*y2))));
    }
} // faddeeva_erf

/******************************************************************************/
// erfi(z) = -i erf(iz)
cmplx faddeeva_erfi(cmplx z, double relerr)
{
    cmplx e = faddeeva_erf(C(-cimag(z),creal(z)), relerr);
    return C(cimag(e), -creal(e));
}

/******************************************************************************/
// erfi(x) = -i erf(ix)
double faddeeva_erfi_re(double x)
{
    return x*x > 720 ? (x > 0 ? Inf : -Inf)
        : exp(x*x) * im_w_of_z(x);
}

/******************************************************************************/
// erfc(x) = 1 - erf(x)
double faddeeva_erfc_re(double x)
{
    return erfc(x); // C99 supplies erfc in math.h
}

/******************************************************************************/
// erfc(z) = 1 - erf(z)
cmplx faddeeva_erfc(cmplx z, double relerr)
{
    double x = creal(z), y = cimag(z);

    if (x == 0.)
        return C(1,
                 /* handle y -> Inf limit manually, since
                    exp(y^2) -> Inf but Im[w(y)] -> 0, so
                    IEEE will give us a NaN when it should be Inf */
                 y*y > 720 ? (y > 0 ? -Inf : Inf)
                 : -exp(y*y) * im_w_of_z(y));
    if (y == 0.) {
        if (x*x > 750) // underflow
            return C(x >= 0 ? 0.0 : 2.0,
                     -y); // preserve sign of 0
        return C(x >= 0 ? exp(-x*x) * faddeeva_erfcx_re(x) 
                 : 2. - exp(-x*x) * faddeeva_erfcx_re(-x),
                 -y); // preserve sign of zero
    }

    double mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
    double mIm_z2 = -2*x*y; // Im(-z^2)
    if (mRe_z2 < -750) // underflow
        return (x >= 0 ? 0.0 : 2.0);

    if (x >= 0)
        return cexp(C(mRe_z2, mIm_z2))
            * faddeeva(C(-y,x), relerr);
    else
        return 2.0 - cexp(C(mRe_z2, mIm_z2))
            * faddeeva(C(y,-x), relerr);
} // faddeeva_erfc

/******************************************************************************/
// compute Dawson(x) = sqrt(pi)/2  *  exp(-x^2) * erfi(x)
double dawson_re(double x)
{
    const double spi2 = 0.8862269254527580136490837416705725913990; // sqrt(pi)/2
    return spi2 * im_w_of_z(x);
} // dawson_re

/******************************************************************************/
// compute Dawson(z) = sqrt(pi)/2  *  exp(-z^2) * erfi(z)
cmplx dawson(cmplx z, double relerr)
{
    const double spi2 = 0.8862269254527580136490837416705725913990; // sqrt(pi)/2
    double x = creal(z), y = cimag(z);

    // handle axes separately for speed & proper handling of x or y = Inf or NaN
    if (y == 0)
        return C(spi2 * im_w_of_z(x),
                 -y); // preserve sign of 0
    if (x == 0) {
        double y2 = y*y;
        if (y2 < 2.5e-5) { // Taylor expansion
            return C(x, // preserve sign of 0
                     y * (1.
                          + y2 * (0.6666666666666666666666666666666666666667
                                  + y2 * 0.26666666666666666666666666666666666667)));
        }
        return C(x, // preserve sign of 0
                 spi2 * (y >= 0 
                         ? exp(y2) - faddeeva_erfcx_re(y)
                         : faddeeva_erfcx_re(-y) - exp(y2)));
    }

    double mRe_z2 = (y - x) * (x + y); // Re(-z^2), being careful of overflow
    double mIm_z2 = -2*x*y; // Im(-z^2)
    cmplx mz2 = C(mRe_z2, mIm_z2); // -z^2

    /* Handle positive and negative x via different formulas,
       using the mirror symmetries of w, to avoid overflow/underflow
       problems from multiplying exponentially large and small quantities. */
    if (y >= 0) {
        if (y < 5e-3) {
            if (fabs(x) < 5e-3)
                goto taylor;
            else if (fabs(mIm_z2) < 5e-3)
                goto taylor_realaxis;
        }
        cmplx res = cexp(mz2) - faddeeva(z, relerr);
        return spi2 * C(-cimag(res), creal(res));
    }
    else { // y < 0
        if (y > -5e-3) { // duplicate from above to avoid fabs(x) call
            if (fabs(x) < 5e-3)
                goto taylor;
            else if (fabs(mIm_z2) < 5e-3)
                goto taylor_realaxis;
        }
        else if (isnan(y))
            return C(x == 0 ? 0 : NaN, NaN);
        cmplx res = faddeeva(-z, relerr) - cexp(mz2);
        return spi2 * C(-cimag(res), creal(res));
    }

    // Use Taylor series for small |z|, to avoid cancellation inaccuracy
    //     dawson(z) = z - 2/3 z^3 + 4/15 z^5 + ...
taylor:
    return z * (1.
                + mz2 * (0.6666666666666666666666666666666666666667
                         + mz2 * 0.2666666666666666666666666666666666666667));

    /* for small |y| and small |xy|, 
       use Taylor series to avoid cancellation inaccuracy:
       dawson(x + iy)
       = D + y^2 (D + x - 2Dx^2)
       + y^4 (D/2 + 5x/6 - 2Dx^2 - x^3/3 + 2Dx^4/3)
       + iy [ (1-2Dx) + 2/3 y^2 (1 - 3Dx - x^2 + 2Dx^3)
       + y^4/15 (4 - 15Dx - 9x^2 + 20Dx^3 + 2x^4 - 4Dx^5) ] + ...
       where D = dawson(x) 

       However, for large |x|, 2Dx -> 1 which gives cancellation problems in
       this series (many of the leading terms cancel).  So, for large |x|,
       we need to substitute a continued-fraction expansion for D.

       dawson(x) = 0.5 / (x-0.5/(x-1/(x-1.5/(x-2/(x-2.5/(x...))))))

       The 6 terms shown here seems to be the minimum needed to be
       accurate as soon as the simpler Taylor expansion above starts
       breaking down.  Using this 6-term expansion, factoring out the
       denominator, and simplifying with Maple, we obtain:

       Re dawson(x + iy) * (-15 + 90x^2 - 60x^4 + 8x^6) / x
       = 33 - 28x^2 + 4x^4 + y^2 (18 - 4x^2) + 4 y^4
       Im dawson(x + iy) * (-15 + 90x^2 - 60x^4 + 8x^6) / y
       = -15 + 24x^2 - 4x^4 + 2/3 y^2 (6x^2 - 15) - 4 y^4

       Finally, for |x| > 5e7, we can use a simpler 1-term continued-fraction
       expansion for the real part, and a 2-term expansion for the imaginary
       part.  (This avoids overflow problems for huge |x|.)  This yields:
     
       Re dawson(x + iy) = [1 + y^2 (1 + y^2/2 - (xy)^2/3)] / (2x)
       Im dawson(x + iy) = y [ -1 - 2/3 y^2 + y^4/15 (2x^2 - 4) ] / (2x^2 - 1)

    */
taylor_realaxis:
    {
        double x2 = x*x;
        if (x2 > 1600) { // |x| > 40
            double y2 = y*y;
            if (x2 > 25e14) {// |x| > 5e7
                double xy2 = (x*y)*(x*y);
                return C((0.5 + y2 * (0.5 + 0.25*y2
                                      - 0.16666666666666666667*xy2)) / x,
                         y * (-1 + y2 * (-0.66666666666666666667
                                         + 0.13333333333333333333*xy2
                                         - 0.26666666666666666667*y2))
                         / (2*x2 - 1));
            }
            return (1. / (-15 + x2*(90 + x2*(-60 + 8*x2)))) *
                C(x * (33 + x2 * (-28 + 4*x2)
                       + y2 * (18 - 4*x2 + 4*y2)),
                  y * (-15 + x2 * (24 - 4*x2)
                       + y2 * (4*x2 - 10 - 4*y2)));
        }
        else {
            double D = spi2 * im_w_of_z(x);
            double y2 = y*y;
            return C
                (D + y2 * (D + x - 2*D*x2)
                 + y2*y2 * (D * (0.5 - x2 * (2 - 0.66666666666666666667*x2))
                            + x * (0.83333333333333333333
                                   - 0.33333333333333333333 * x2)),
                 y * (1 - 2*D*x
                      + y2 * 0.66666666666666666667 * (1 - x2 - D*x * (3 - 2*x2))
                      + y2*y2 * (0.26666666666666666667 -
                                 x2 * (0.6 - 0.13333333333333333333 * x2)
                                 - D*x * (1 - x2 * (1.3333333333333333333
                                                    - 0.26666666666666666667 * x2)))));
        }
    }
} // dawson

/******************************************************************************/
// return sinc(x) = sin(x)/x, given both x and sin(x) 
// [since we only use this in cases where sin(x) has already been computed]
static inline double sinc(double x, double sinx) { 
    return fabs(x) < 1e-4 ? 1 - (0.1666666666666666666667)*x*x : sinx / x; 
}

/******************************************************************************/
// sinh(x) via Taylor series, accurate to machine precision for |x| < 1e-2
static inline double sinh_taylor(double x) {
    return x * (1 + (x*x) * (0.1666666666666666666667
                             + 0.00833333333333333333333 * (x*x)));
}

/******************************************************************************/
static inline double sqr(double x) { return x*x; }

/******************************************************************************/
// precomputed table of expa2n2[n-1] = exp(-a2*n*n)
// for double-precision a2 = 0.26865... in faddeeva, below.
static const double expa2n2[] = {
    7.64405281671221563e-01,
    3.41424527166548425e-01,
    8.91072646929412548e-02,
    1.35887299055460086e-02,
    1.21085455253437481e-03,
    6.30452613933449404e-05,
    1.91805156577114683e-06,
    3.40969447714832381e-08,
    3.54175089099469393e-10,
    2.14965079583260682e-12,
    7.62368911833724354e-15,
    1.57982797110681093e-17,
    1.91294189103582677e-20,
    1.35344656764205340e-23,
    5.59535712428588720e-27,
    1.35164257972401769e-30,
    1.90784582843501167e-34,
    1.57351920291442930e-38,
    7.58312432328032845e-43,
    2.13536275438697082e-47,
    3.51352063787195769e-52,
    3.37800830266396920e-57,
    1.89769439468301000e-62,
    6.22929926072668851e-68,
    1.19481172006938722e-73,
    1.33908181133005953e-79,
    8.76924303483223939e-86,
    3.35555576166254986e-92,
    7.50264110688173024e-99,
    9.80192200745410268e-106,
    7.48265412822268959e-113,
    3.33770122566809425e-120,
    8.69934598159861140e-128,
    1.32486951484088852e-135,
    1.17898144201315253e-143,
    6.13039120236180012e-152,
    1.86258785950822098e-160,
    3.30668408201432783e-169,
    3.43017280887946235e-178,
    2.07915397775808219e-187,
    7.36384545323984966e-197,
    1.52394760394085741e-206,
    1.84281935046532100e-216,
    1.30209553802992923e-226,
    5.37588903521080531e-237,
    1.29689584599763145e-247,
    1.82813078022866562e-258,
    1.50576355348684241e-269,
    7.24692320799294194e-281,
    2.03797051314726829e-292,
    3.34880215927873807e-304,
    0.0 // underflow (also prevents reads past array end, below)
}; // expa2n2

/******************************************************************************/
cmplx faddeeva(cmplx z, double relerr)
{
    if (creal(z) == 0.0)
        return C(faddeeva_erfcx_re(cimag(z)), creal(z));
               // give correct sign of 0 in cimag(w)
    else if (cimag(z) == 0)
        return C(exp(-sqr(creal(z))),  im_w_of_z(creal(z)));

    double a, a2, c;
    if (relerr <= DBL_EPSILON) {
        relerr = DBL_EPSILON;
        a = 0.518321480430085929872; // pi / sqrt(-log(eps*0.5))
        c = 0.329973702884629072537; // (2/pi) * a;
        a2 = 0.268657157075235951582; // a^2
    }
    else {
        const double pi = 3.14159265358979323846264338327950288419716939937510582;
        if (relerr > 0.1) relerr = 0.1; // not sensible to compute < 1 digit
        a = pi / sqrt(-log(relerr*0.5));
        c = (2/pi)*a;
        a2 = a*a;
    }
    const double x = fabs(creal(z));
    const double y = cimag(z), ya = fabs(y);

    cmplx ret = 0.; // return value

    double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0;

#define USE_CONTINUED_FRACTION 1 // 1 to use continued fraction for large |z|

#if USE_CONTINUED_FRACTION
    if (ya > 7 || (x > 6  // continued fraction is faster
                   /* As pointed out by M. Zaghloul, the continued
                      fraction seems to give a large relative error in
                      Re w(z) for |x| ~ 6 and small |y|, so use
                      algorithm 816 in this region: */
                   && (ya > 0.1 || (x > 8 && ya > 1e-10) || x > 28))) {
    
        /* Poppe & Wijers suggest using a number of terms
           nu = 3 + 1442 / (26*rho + 77)
           where rho = sqrt((x/x0)^2 + (y/y0)^2) where x0=6.3, y0=4.4.
           (They only use this expansion for rho >= 1, but rho a little less
           than 1 seems okay too.)
           Instead, I did my own fit to a slightly different function
           that avoids the hypotenuse calculation, using NLopt to minimize
           the sum of the squares of the errors in nu with the constraint
           that the estimated nu be >= minimum nu to attain machine precision.
           I also separate the regions where nu == 2 and nu == 1. */
        const double ispi = 0.56418958354775628694807945156; // 1 / sqrt(pi)
        double xs = y < 0 ? -creal(z) : creal(z); // compute for -z if y < 0
        if (x + ya > 4000) { // nu <= 2
            if (x + ya > 1e7) { // nu == 1, w(z) = i/sqrt(pi) / z
                // scale to avoid overflow
                if (x > ya) {
                    double yax = ya / xs; 
                    double denom = ispi / (xs + yax*ya);
                    ret = C(denom*yax, denom);
                }
                else if (isinf(ya))
                    return ((isnan(x) || y < 0) 
                            ? C(NaN,NaN) : C(0,0));
                else {
                    double xya = xs / ya;
                    double denom = ispi / (xya*xs + ya);
                    ret = C(denom, denom*xya);
                }
            }
            else { // nu == 2, w(z) = i/sqrt(pi) * z / (z*z - 0.5)
                double dr = xs*xs - ya*ya - 0.5, di = 2*xs*ya;
                double denom = ispi / (dr*dr + di*di);
                ret = C(denom * (xs*di-ya*dr), denom * (xs*dr+ya*di));
            }
        }
        else { // compute nu(z) estimate and do general continued fraction
            const double c0=3.9, c1=11.398, c2=0.08254, c3=0.1421, c4=0.2023; // fit
            double nu = floor(c0 + c1 / (c2*x + c3*ya + c4));
            double wr = xs, wi = ya;
            for (nu = 0.5 * (nu - 1); nu > 0.4; nu -= 0.5) {
                // w <- z - nu/w:
                double denom = nu / (wr*wr + wi*wi);
                wr = xs - wr * denom;
                wi = ya + wi * denom;
            }
            { // w(z) = i/sqrt(pi) / w:
                double denom = ispi / (wr*wr + wi*wi);
                ret = C(denom*wi, denom*wr);
            }
        }
        if (y < 0) {
            // use w(z) = 2.0*exp(-z*z) - w(-z), 
            // but be careful of overflow in exp(-z*z) 
            //                                = exp(-(xs*xs-ya*ya) -2*i*xs*ya) 
            return 2.0*cexp(C((ya-xs)*(xs+ya), 2*xs*y)) - ret;
        }
        else
            return ret;
    }
#else // !USE_CONTINUED_FRACTION
    if (x + ya > 1e7) { // w(z) = i/sqrt(pi) / z, to machine precision
        const double ispi = 0.56418958354775628694807945156; // 1 / sqrt(pi)
        double xs = y < 0 ? -creal(z) : creal(z); // compute for -z if y < 0
        // scale to avoid overflow
        if (x > ya) {
            double yax = ya / xs; 
            double denom = ispi / (xs + yax*ya);
            ret = C(denom*yax, denom);
        }
        else {
            double xya = xs / ya;
            double denom = ispi / (xya*xs + ya);
            ret = C(denom, denom*xya);
        }
        if (y < 0) {
            // use w(z) = 2.0*exp(-z*z) - w(-z), 
            // but be careful of overflow in exp(-z*z) 
            //                                = exp(-(xs*xs-ya*ya) -2*i*xs*ya) 
            return 2.0*cexp(C((ya-xs)*(xs+ya), 2*xs*y)) - ret;
        }
        else
            return ret;
    }
#endif // !USE_CONTINUED_FRACTION 

    /* Note: The test that seems to be suggested in the paper is x <
       sqrt(-log(DBL_MIN)), about 26.6, since otherwise exp(-x^2)
       underflows to zero and sum1,sum2,sum4 are zero.  However, long
       before this occurs, the sum1,sum2,sum4 contributions are
       negligible in double precision; I find that this happens for x >
       about 6, for all y.  On the other hand, I find that the case
       where we compute all of the sums is faster (at least with the
       precomputed expa2n2 table) until about x=10.  Furthermore, if we
       try to compute all of the sums for x > 20, I find that we
       sometimes run into numerical problems because underflow/overflow
       problems start to appear in the various coefficients of the sums,
       below.  Therefore, we use x < 10 here. */
    else if (x < 10) {
        double prod2ax = 1, prodm2ax = 1;
        double expx2;

        if (isnan(y))
            return C(y,y);
    
        /* Somewhat ugly copy-and-paste duplication here, but I see significant
           speedups from using the special-case code with the precomputed
           exponential, and the x < 5e-4 special case is needed for accuracy. */

        if (relerr == DBL_EPSILON) { // use precomputed exp(-a2*(n*n)) table
            if (x < 5e-4) { // compute sum4 and sum5 together as sum5-sum4
                const double x2 = x*x;
                expx2 = 1 - x2 * (1 - 0.5*x2); // exp(-x*x) via Taylor
                // compute exp(2*a*x) and exp(-2*a*x) via Taylor, to double precision
                const double ax2 = 1.036642960860171859744*x; // 2*a*x
                const double exp2ax =
                    1 + ax2 * (1 + ax2 * (0.5 + 0.166666666666666666667*ax2));
                const double expm2ax =
                    1 - ax2 * (1 - ax2 * (0.5 - 0.166666666666666666667*ax2));
                for (int n = 1; 1; ++n) {
                    const double coef = expa2n2[n-1] * expx2 / (a2*(n*n) + y*y);
                    prod2ax *= exp2ax;
                    prodm2ax *= expm2ax;
                    sum1 += coef;
                    sum2 += coef * prodm2ax;
                    sum3 += coef * prod2ax;
          
                    // really = sum5 - sum4
                    sum5 += coef * (2*a) * n * sinh_taylor((2*a)*n*x);
          
                    // test convergence via sum3
                    if (coef * prod2ax < relerr * sum3) break;
                }
            }
            else { // x > 5e-4, compute sum4 and sum5 separately
                expx2 = exp(-x*x);
                const double exp2ax = exp((2*a)*x), expm2ax = 1 / exp2ax;
                for (int n = 1; 1; ++n) {
                    const double coef = expa2n2[n-1] * expx2 / (a2*(n*n) + y*y);
                    prod2ax *= exp2ax;
                    prodm2ax *= expm2ax;
                    sum1 += coef;
                    sum2 += coef * prodm2ax;
                    sum4 += (coef * prodm2ax) * (a*n);
                    sum3 += coef * prod2ax;
                    sum5 += (coef * prod2ax) * (a*n);
                    // test convergence via sum5, since this sum has the slowest decay
                    if ((coef * prod2ax) * (a*n) < relerr * sum5) break;
                }
            }
        }
        else { // relerr != DBL_EPSILON, compute exp(-a2*(n*n)) on the fly
            const double exp2ax = exp((2*a)*x), expm2ax = 1 / exp2ax;
            if (x < 5e-4) { // compute sum4 and sum5 together as sum5-sum4
                const double x2 = x*x;
                expx2 = 1 - x2 * (1 - 0.5*x2); // exp(-x*x) via Taylor
                for (int n = 1; 1; ++n) {
                    const double coef = exp(-a2*(n*n)) * expx2 / (a2*(n*n) + y*y);
                    prod2ax *= exp2ax;
                    prodm2ax *= expm2ax;
                    sum1 += coef;
                    sum2 += coef * prodm2ax;
                    sum3 += coef * prod2ax;
          
                    // really = sum5 - sum4
                    sum5 += coef * (2*a) * n * sinh_taylor((2*a)*n*x);
          
                    // test convergence via sum3
                    if (coef * prod2ax < relerr * sum3) break;
                }
            }
            else { // x > 5e-4, compute sum4 and sum5 separately
                expx2 = exp(-x*x);
                for (int n = 1; 1; ++n) {
                    const double coef = exp(-a2*(n*n)) * expx2 / (a2*(n*n) + y*y);
                    prod2ax *= exp2ax;
                    prodm2ax *= expm2ax;
                    sum1 += coef;
                    sum2 += coef * prodm2ax;
                    sum4 += (coef * prodm2ax) * (a*n);
                    sum3 += coef * prod2ax;
                    sum5 += (coef * prod2ax) * (a*n);
                    // test convergence via sum5, since this sum has the slowest decay
                    if ((coef * prod2ax) * (a*n) < relerr * sum5) break;
                }
            }
        }
        const double expx2erfcxy = // avoid spurious overflow for large negative y
            y > -6 // for y < -6, erfcx(y) = 2*exp(y*y) to double precision
            ? expx2*faddeeva_erfcx_re(y) : 2*exp(y*y-x*x);
        if (y > 5) { // imaginary terms cancel
            const double sinxy = sin(x*y);
            ret = (expx2erfcxy - c*y*sum1) * cos(2*x*y)
                + (c*x*expx2) * sinxy * sinc(x*y, sinxy);
        }
        else {
            double xs = creal(z);
            const double sinxy = sin(xs*y);
            const double sin2xy = sin(2*xs*y), cos2xy = cos(2*xs*y);
            const double coef1 = expx2erfcxy - c*y*sum1;
            const double coef2 = c*xs*expx2;
            ret = C(coef1 * cos2xy + coef2 * sinxy * sinc(xs*y, sinxy),
                    coef2 * sinc(2*xs*y, sin2xy) - coef1 * sin2xy);
        }
    }
    else { // x large: only sum3 & sum5 contribute (see above note)    
        if (isnan(x))
            return C(x,x);
        if (isnan(y))
            return C(y,y);

#if USE_CONTINUED_FRACTION
        ret = exp(-x*x); // |y| < 1e-10, so we only need exp(-x*x) term
#else
        if (y < 0) {
            /* erfcx(y) ~ 2*exp(y*y) + (< 1) if y < 0, so
               erfcx(y)*exp(-x*x) ~ 2*exp(y*y-x*x) term may not be negligible
               if y*y - x*x > -36 or so.  So, compute this term just in case.
               We also need the -exp(-x*x) term to compute Re[w] accurately
               in the case where y is very small. */
            ret = cpolar(2*exp(y*y-x*x) - exp(-x*x), -2*creal(z)*y);
        }
        else
            ret = exp(-x*x); // not negligible in real part if y very small
#endif
        // (round instead of ceil as in original paper; note that x/a > 1 here)
        double n0 = floor(x/a + 0.5); // sum in both directions, starting at n0
        double dx = a*n0 - x;
        sum3 = exp(-dx*dx) / (a2*(n0*n0) + y*y);
        sum5 = a*n0 * sum3;
        double exp1 = exp(4*a*dx), exp1dn = 1;
        int dn;
        for (dn = 1; n0 - dn > 0; ++dn) { // loop over n0-dn and n0+dn terms
            double np = n0 + dn, nm = n0 - dn;
            double tp = exp(-sqr(a*dn+dx));
            double tm = tp * (exp1dn *= exp1); // trick to get tm from tp
            tp /= (a2*(np*np) + y*y);
            tm /= (a2*(nm*nm) + y*y);
            sum3 += tp + tm;
            sum5 += a * (np * tp + nm * tm);
            if (a * (np * tp + nm * tm) < relerr * sum5) goto finish;
        }
        while (1) { // loop over n0+dn terms only (since n0-dn <= 0)
            double np = n0 + dn++;
            double tp = exp(-sqr(a*dn+dx)) / (a2*(np*np) + y*y);
            sum3 += tp;
            sum5 += a * np * tp;
            if (a * np * tp < relerr * sum5) goto finish;
        }
    }
finish:
    return ret + C((0.5*c)*y*(sum2+sum3), 
                   (0.5*c)*copysign(sum5-sum4, creal(z)));
} // faddeeva
