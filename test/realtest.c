/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File realtest.c
 *   Compare real function and their complex counterparts for various real arguments
 *
 * Copyright:
 *   (C) 2012 Massachusetts Institute of Technology
 *   (C) 2013 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   ../LICENSE
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
 */

#include "cerf.h"
#include "testtool.h"

// For testing a real function against its complex counterpart.

void real_tests(
    result_t* result, const char* name, _cerf_cmplx (*F)(_cerf_cmplx), double (*FRE)(double),
    double xmin, double xmax)
{
    char info[30];

    // Arguments +-x from logarithmic grid.
    const int n=10000;
    for (int i = 0; i < n; ++i) {
        double x = pow(10., -300. + i * 600. / (n - 1));
        if (x<xmin || x>xmax)
            continue;
        snprintf(info, 30, "%s(%g)", name, x);
        rtest(result, 1e-13, creal(F(C(x, 0))), FRE(x), info);
        rtest(result, 1e-13, creal(F(C(-x, 0))), FRE(-x), info);
	// Also test continuity if complex argument has a tiny imaginary part:
        rtest(result, 1e-10, creal(F(C(x, x * 1e-10))), FRE(x), info);
        rtest(result, 1e-10, creal(F(C(-x, x * 1e-10))), FRE(-x), info);
        rtest(result, 1e-6, creal(F(C(x, x * 1e-6))), FRE(x), info);
        rtest(result, 1e-6, creal(F(C(-x, x * 1e-6))), FRE(-x), info);
    }

    // Arguments 0, +-inf, nan.
    snprintf(info, 30, "0");
    rtest(result, 1e-13, creal(F(C(0., 0.))), FRE(0.), info );
    snprintf(info, 30, "%s(Inf)", name);
    rtest(result, 1e-13, creal(F(C(Inf, 0.))), FRE(Inf), info );
    snprintf(info, 30, "%s(-Inf)", name);
    rtest(result, 1e-13, creal(F(C(-Inf, 0.))), FRE(-Inf), info);
    snprintf(info, 30, "%s(NaN)", name);
    rtest(result, 1e-13, creal(F(C(NaN, 0.))), FRE(NaN), info);
}

int main(void)
{
    result_t result = {0, 0};

    real_tests(&result, "erf", cerf, erf, 1e-300, 1e300);
    real_tests(&result, "erfi", cerfi, erfi, 1e-300, 1e300);
    real_tests(&result, "erfc", cerfc, erfc, 1e-300, 1e300);
    real_tests(&result, "erfcx", cerfcx, erfcx, 1e-300, 1e300);
    real_tests(&result, "dawson", cdawson, dawson, 1e-300, 1e150);

    printf("%i/%i tests failed\n", result.failed, result.total);
    return result.failed;
}
