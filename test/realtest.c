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

const double errBound = 1e-13;

// For testing the Dawson and error functions for the special case of a real argument

void xTest(
    result_t* result, const char* name, _cerf_cmplx (*F)(_cerf_cmplx), double (*FRE)(double),
    double isc)
{
    char info[30];
    for (int i = 0; i < 10000; ++i) {
        double x = pow(10., -300. + i * 600. / (10000 - 1));
        snprintf(info, 30, "%s(%g)", name, x);
        rtest(result, 1e-13, creal(F(C(x, x * isc))), FRE(x), info);
    }
}

// For testing the Dawson and error functions for the special case of an infinite argument

void iTest(result_t* result, const char* name, _cerf_cmplx (*F)(_cerf_cmplx), double (*FRE)(double))
{
    rtest(result, 1e-13, FRE(Inf), creal(F(C(Inf, 0.))), name );
    rtest(result, 1e-13, FRE(-Inf), creal(F(C(-Inf, 0.))), name);
    rtest(result, 1e-13, FRE(NaN), creal(F(C(NaN, 0.))), name);
}

int main(void)
{
    result_t result = {0, 0};

    xTest(&result, "erf", cerf, erf, 1e-20);
    iTest(&result, "erf", cerf, erf);

    xTest(&result, "erfi", cerfi, erfi, 0);
    iTest(&result, "erfi", cerfi, erfi);

    xTest(&result, "erfc", cerfc, erfc, 1e-20);
    iTest(&result, "erfc", cerfc, erfc);

    xTest(&result, "erfcx", cerfcx, erfcx, 0);
    iTest(&result, "erfcx", cerfcx, erfcx);

    xTest(&result, "dawson", cdawson, dawson, 1e-20);
    iTest(&result, "dawson", cdawson, dawson);

    printf("%i/%i tests failed", result.failed, result.total);
    return result.failed;
}
