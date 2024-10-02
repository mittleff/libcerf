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
 *   Steven G. Johnson, Massachusetts Institute of Technology, 2012, core author
 *   Joachim Wuttke, Forschungszentrum Jülich, 2013, package maintainer
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

static inline double sinh_taylor(double x) {
    // sinh(x) via Taylor series, accurate to machine precision for |x| < 1e-2
    return x * (1 + (x * x) * (0.1666666666666666666667 +
                               0.00833333333333333333333 * (x * x)));
}

static inline double sqr(double x) { return x * x; }

/******************************************************************************/
/* precomputed table of expa2n2[n-1] = exp(-a2*n*n)                           */
/* for double-precision a2 = 0.26865... in w_of_z, below.                     */
/******************************************************************************/

static const double expa2n2[] = {
    7.6440528167122157e-1,
    3.41424527166548419e-1,
    8.91072646929412376e-2,
    1.35887299055460053e-2,
    1.21085455253437473e-3,
    6.30452613933448798e-5,
    1.91805156577114627e-6,
    3.40969447714832129e-8,
    3.54175089099468534e-10,
    2.14965079583260701e-12,
    7.62368911833724214e-15,
    1.57982797110680523e-17,
    1.91294189103582847e-20,
    1.35344656764205201e-23,
    5.59535712428587329e-27,
    1.35164257972401336e-30,
    1.90784582843499203e-34,
    1.573519202914414e-38,
    7.58312432328031747e-43,
    2.13536275438697177e-47,
    3.51352063787194301e-52,
    3.37800830266396575e-57,
    1.89769439468300171e-62,
    6.22929926072660027e-68,
    1.19481172006938479e-73,
    1.33908181133006436e-79,
    8.76924303483226468e-86,
    3.35555576166253504e-92,
    7.5026411068815959e-99,
    9.80192200745400666e-106,
    7.48265412822263025e-113,
    3.33770122566805208e-120,
    8.69934598159840512e-128,
    1.3248695148408338e-135,
    1.17898144201314251e-143,
    6.13039120236156112e-152,
    1.86258785950818541e-160,
    3.30668408201430881e-169,
    3.43017280887946632e-178,
    2.07915397775808552e-187,
    7.36384545323981754e-197,
    1.52394760394083166e-206,
    1.84281935046525516e-216,
    1.30209553802992364e-226,
    5.37588903521091667e-237,
    1.29689584599760859e-247,
    1.82813078022865549e-258,
    1.50576355348675694e-269,
    7.24692320799252486e-281,
    2.03797051314725175e-292,
    3.34880215927866455e-304,
    0.0 // underflow (will force termination of loops, and thereby prevent reading past array end)
}; // expa2n2

/******************************************************************************/
/*  w_of_z, Faddeeva's scaled complex error function                          */
/******************************************************************************/

_cerf_cmplx w_of_z(_cerf_cmplx z) {
#ifdef CERF_INTROSPECT
    cerf_algorithm = -1;
    cerf_nofterms = 0;
#endif

    // Steven G. Johnson, October 2012.

    if (creal(z) == 0.0) {
#ifdef CERF_INTROSPECT
        cerf_algorithm = 400;
#endif
        // Purely imaginary input, purely real output.
        // However, use creal(z) to give correct sign of 0 in cimag(w).
        return C(erfcx(cimag(z)), creal(z));
    }
    if (cimag(z) == 0) {
#ifdef CERF_INTROSPECT
        cerf_algorithm = 500;
#endif
        // Purely real input, complex output.
        // Avoid floating underflow for real term of large z.
        const double Wreal = fabs(creal(z)) > 27. ? 0. : exp(-sqr(creal(z)));
        const double Wimag = im_w_of_x(creal(z));
        return C(Wreal, Wimag);
    }

    const double relerr = DBL_EPSILON;
    const double a = 0.518321480430085929872;  // pi / sqrt(-log(eps*0.5))
    const double c = 0.329973702884629072537;  // (2/pi) * a;
    const double a2 = 0.268657157075235951582; // a^2

    const double x = fabs(creal(z));
    const double y = cimag(z);
    const double ya = fabs(y);

    _cerf_cmplx ret = 0.; // return value

    double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0;

    if (ya > 7 || (x > 6 && (ya > 0.1 || (x > 8 && ya > 1e-10) || x > 28))) {
        /*
	   Continued-fraction expansion,
           similar to those described by Gautschi (1970) and Poppe & Wijers (1990).

           Prefered for large z because it is fast.

           However, as pointed out by M. Zaghloul, the continued fraction seems to
           give a large relative error in Re w(z) for |x| ~ 6 and small |y|.
           In this region, excluded by the additional conditions in the if clause above,
           we fall back to ACM algorithm 916 below.

           Poppe & Wijers suggest using a number of terms
           nu = 3 + 1442 / (26*rho + 77)
           where rho = sqrt((x/x0)^2 + (y/y0)^2) where x0=6.3, y0=4.4.
           (They only use this expansion for rho >= 1, but rho a little less
           than 1 seems okay too.)
           Instead, I [SGJ] did my own fit to a slightly different function
           that avoids the hypotenuse calculation, using NLopt to minimize
           the sum of the squares of the errors in nu with the constraint
           that the estimated nu be >= minimum nu to attain machine precision.
           I also separate the regions where nu == 2 and nu == 1.
	*/

        const double ispi = 0.56418958354775628694807945156; // 1 / sqrt(pi)
        const double xs = y < 0 ? -creal(z) : creal(z); // compute for -z if y < 0

        if (x + ya > 4000) {                      // nu <= 2
            if (x + ya > 1e7) {                     // nu == 1, w(z) = i/sqrt(pi) / z
                // scale to avoid overflow
                if (x > ya) {
#ifdef CERF_INTROSPECT
                    cerf_algorithm = 101;
#endif
                    const double yax = ya / xs;
                    const double denom = ispi / (xs + yax*ya);
                    ret = C(denom*yax, denom);
                } else if (isinf(ya)) {
#ifdef CERF_INTROSPECT
                    cerf_algorithm = 102;
#endif
                    return ((isnan(x) || y < 0) ? C(NaN, NaN) : C(0, 0));
                } else {
#ifdef CERF_INTROSPECT
                    cerf_algorithm = 103;
#endif
                    const double xya = xs / ya;
                    const double denom = ispi / (xya*xs + ya);
                    ret = C(denom, denom*xya);
                }

            } else { // nu == 2, w(z) = i/sqrt(pi)*z / (z*z - 0.5)
#ifdef CERF_INTROSPECT
                cerf_algorithm = 104;
#endif
                const double dr = xs*xs - ya*ya - 0.5, di = 2*xs*ya;
                const double denom = ispi / (dr*dr + di*di);
                ret = C(denom * (xs*di - ya*dr), denom * (xs*dr + ya*di));
            }

        } else {
            // compute nu(z) estimate and do general continued fraction
#ifdef CERF_INTROSPECT
            cerf_algorithm = 105;
#endif
            const double c0 = 3.9, c1 = 11.398, c2 = 0.08254, c3 = 0.1421,
                c4 = 0.2023; // fit
            double nu = floor(c0 + c1 / (c2*x + c3*ya + c4));
            double wr = xs;
            double wi = ya;
            for (nu = 0.5 * (nu - 1); nu > 0.4; nu -= 0.5) {
                // w <- z - nu/w:
                double denom = nu / (wr*wr + wi*wi);
                wr = xs - wr*denom;
                wi = ya + wi*denom;
            }
            // w(z) = i/sqrt(pi) / w:
            double denom = ispi / (wr*wr + wi*wi);
            ret = C(denom*wi, denom*wr);
        }

        if (y < 0) {
            /*
	       use w(z) = 2.0*exp(-z*z) - w(-z),
               but be careful of overflow in exp(-z*z) = exp(-(xs*xs-ya*ya) -2*i*xs*ya)
	    */
#ifdef CERF_INTROSPECT
            cerf_algorithm += 10;
#endif
            return 2.0 * cexp(C((ya - xs) * (xs + ya), 2*xs*y)) - ret;
        } else
            return ret;
    }

    else if (x < 10) {
	/*
	   ACM algorithm 916 by Zaghloul & Ali (2011), which is generally competitive
           at small |z|, and more accurate than the Poppe & Wijers expansion in some
           regions, e.g. in the vicinity of z=1+i.

       	   Applicability range: The paper seems to suggest x < sqrt(-log(DBL_MIN)),
           about 26.6, since otherwise exp(-x^2) underflows to zero and sum1,sum2,sum4
           are zero.  However, long before this occurs, the sum1,sum2,sum4 contributions
           are negligible in double precision; I [SGJ] find that this happens for
           x > about 6, for all y.  On the other hand, I find that the case
       	   where we compute all of the sums is faster (at least with the
       	   precomputed expa2n2 table) until about x=10.
        */

	double prod2ax = 1, prodm2ax = 1;
	double expx2;

        if (isnan(y)) {
#ifdef CERF_INTROSPECT
            cerf_algorithm = 299;
#endif
            return C(y, y);
        }

        if (x < 5e-4) {
	    /*
	       compute sum4 and sum5 together as sum5-sum4
               This special case is needed for accuracy.
	    */
#ifdef CERF_INTROSPECT
            cerf_algorithm = 201;
#endif
            const double x2 = x*x;
            expx2 = 1 - x2 * (1 - 0.5*x2); // exp(-x*x) via Taylor
            // compute exp(2*a*x) and exp(-2*a*x) via Taylor, to double precision
            const double ax2 = 1.036642960860171859744 * x; // 2*a*x
            const double exp2ax =
                1 + ax2 * (1 + ax2 * (0.5 + 0.166666666666666666667 * ax2));
            const double expm2ax =
                1 - ax2 * (1 - ax2 * (0.5 - 0.166666666666666666667 * ax2));
            for (int n = 1;; ++n) {
                const double coef = expa2n2[n - 1] * expx2 / (a2 * (n*n) + y*y);
                prod2ax *= exp2ax;
                prodm2ax *= expm2ax;
                sum1 += coef;
                sum2 += coef * prodm2ax;
                sum3 += coef * prod2ax;

                // really = sum5 - sum4
                sum5 += coef * (2*a) * n * sinh_taylor((2*a) * n * x);

                // test convergence via sum3; for termination rely on coef[n_max] = 0
                if (coef * prod2ax < relerr * sum3) {
#ifdef CERF_INTROSPECT
                    cerf_nofterms = n;
#endif
                    break;
                }
            }

        } else {
            /*
	       x > 5e-4, compute sum4 and sum5 separately
	    */

#ifdef CERF_INTROSPECT
            cerf_algorithm = 202;
#endif
            expx2 = exp(-x*x);
            const double exp2ax = exp((2*a) * x), expm2ax = 1 / exp2ax;
            for (int n = 1;; ++n) {
                const double coef = expa2n2[n - 1] * expx2 / (a2 * (n*n) + y*y);
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
#ifdef CERF_INTROSPECT
                    cerf_nofterms = n;
#endif
                    break;
                }
            }
        }
        const double expx2erfcxy = y < -6  ? 2*exp(y*y-x*x) : expx2*erfcx(y);
           /*
	      The second case has the exact expression.
              In the first case, to avoid spurious overflow for large negative y,
              we approximate erfcx(y) by 2*exp(y^2), which is accurate to double precision.
	   */
	   // TODO: check exact location of cross-over

        if (y > 5) { // imaginary terms cancel
#ifdef CERF_INTROSPECT
            cerf_algorithm += 10;
#endif
            const double sinxy = sin(x*y);
            ret = (expx2erfcxy - c*y*sum1) * cos(2*x*y) +
                (c*x*expx2) * sinxy * sinc(x*y, sinxy);
        } else {
#ifdef CERF_INTROSPECT
            cerf_algorithm += 20;
#endif
            double xs = creal(z);
            const double sinxy = sin(xs*y);
            const double sin2xy = sin(2*xs*y), cos2xy = cos(2*xs*y);
            const double coef1 = expx2erfcxy - c*y*sum1;
            const double coef2 = c*xs*expx2;
            ret = C(coef1*cos2xy + coef2*sinxy*sinc(xs*y, sinxy),
                    coef2*sinc(2*xs*y, sin2xy) - coef1*sin2xy);
        }

    } else {
	/*
	   Currently, this case cannot be reached.

           Still ACM algorithm 916 by Zaghloul & Ali (2011), modified for large x.

	   In the original algorithm, if we try to compute all of the sums for x > 20,
           I [SGJ] find that we sometimes run into numerical problems because
           underflow/overflow problems start to appear in the coefficients of some sums.

           Here, only sum3 & sum5 contribute.
	*/
#ifdef CERF_INTROSPECT
        cerf_algorithm = 300;
#endif

        if (isnan(x))
            return C(x, x);
        if (isnan(y))
            return C(y, y);

        ret = exp(-x*x); // |y| < 1e-10, so we only need exp(-x*x) term
        // (round instead of ceil as in original paper; note that x/a > 1 here)
        const double n0 = floor(x / a + 0.5); // sum in both directions, starting at n0
        const double dx = a*n0 - x;
        sum3 = exp(-dx*dx) / (a2 * (n0*n0) + y*y);
        sum5 = a*n0*sum3;
        const double exp1 = exp(4*a*dx);
        double exp1dn = 1;
        int dn;
        for (dn = 1; n0 - dn > 0; ++dn) { // loop over n0-dn and n0+dn terms
            const double np = n0 + dn, nm = n0 - dn;
            double tp = exp(-sqr(a*dn + dx));
            double tm = tp * (exp1dn *= exp1); // trick to get tm from tp
            tp /= (a2 * (np*np) + y*y);
            tm /= (a2 * (nm*nm) + y*y);
            sum3 += tp + tm;
            sum5 += a * (np*tp + nm*tm);
            if (a * (np*tp + nm*tm) < relerr * sum5)
                goto finish;
        }
        while (1) { // loop over n0+dn terms only (since n0-dn <= 0)
            const double np = n0 + dn++;
            const double tp = exp(-sqr(a*dn + dx)) / (a2 * (np*np) + y*y);
            sum3 += tp;
            sum5 += a*np*tp;
            if (a*np*tp < relerr * sum5)
                goto finish;
        }
    }
finish:
    return ret + C((c/2) * y * (sum2 + sum3), (c/2) * copysign(sum5 - sum4, creal(z)));
} // w_of_z
