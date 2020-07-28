/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File realtest.c
 *   Compare complex and real function for various real arguments
 *
 * Copyright:
 *   (C) 2012 Massachusetts Institute of Technology
 *   (C) 2013 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   ../LICENSE
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
 */

#include "cerf.h"
#include "testtool.h"

#include <float.h>
#include <math.h>
#include <stdio.h>

const double errBound = 1e-13;

/******************************************************************************/
/*  Auxiliary routines                                                        */
/******************************************************************************/

// For testing the Dawson and error functions for the special case of a real argument

void xTest(
    int* fail, const char* fctName, _cerf_cmplx (*F)(_cerf_cmplx), double (*FRE)(double),
    double isc)
{
    printf("############# %s(x) tests #############\n", fctName);
    double errmax = 0;
    for (int i = 0; i < 10000; ++i) {
        double x = pow(10., -300. + i * 600. / (10000 - 1));
        double re_err = relerr(FRE(x), creal(F(C(x, x * isc))));
        if (re_err > errmax)
            errmax = re_err;
        re_err = relerr(FRE(-x), creal(F(C(-x, x * isc))));
        if (re_err > errmax)
            errmax = re_err;
    }
    if (errmax > errBound) {
        printf("FAILURE -- relative error %g too large!\n", errmax);
        ++(*fail);
    } else
        printf("SUCCESS (max relative error = %g)\n", errmax);
}

// For testing the Dawson and error functions for the special case of an infinite argument

void iTest(int* fail, const char* fctName, _cerf_cmplx (*F)(_cerf_cmplx), double (*FRE)(double))
{
    printf("############# %s(inf) tests ###########\n", fctName);
    double errmax = 0;
    double re_err = relerr(FRE(Inf), creal(F(C(Inf, 0.))));
    if (re_err > errmax)
        errmax = re_err;
    re_err = relerr(FRE(-Inf), creal(F(C(-Inf, 0.))));
    if (re_err > errmax)
        errmax = re_err;
    re_err = relerr(FRE(NaN), creal(F(C(NaN, 0.))));
    if (re_err > errmax)
        errmax = re_err;
    if (errmax > errBound) {
        printf("FAILURE -- relative error %g too large!\n", errmax);
        ++(*fail);
    } else
        printf("SUCCESS (max relative error = %g)\n", errmax);
}

/******************************************************************************/
/*  Main: test sequence                                                       */
/******************************************************************************/

int main(void)
{
    int fail = 0;

    xTest(&fail, "erf", cerf, erf, 1e-20);
    iTest(&fail, "erf", cerf, erf);

    xTest(&fail, "erfi", cerfi, erfi, 0);
    iTest(&fail, "erfi", cerfi, erfi);

    xTest(&fail, "erfc", cerfc, erfc, 1e-20);
    iTest(&fail, "erfc", cerfc, erfc);

    xTest(&fail, "erfcx", cerfcx, erfcx, 0);
    iTest(&fail, "erfcx", cerfcx, erfcx);

    xTest(&fail, "dawson", cdawson, dawson, 1e-20);
    iTest(&fail, "dawson", cdawson, dawson);

    printf("#####################################\n");
    if (fail) {
        printf("IN TOTAL, FAILURE IN %i TESTS\n", fail);
        return 1;
    } else {
        printf("OVERALL SUCCESS\n");
        return 0;
    }
}
