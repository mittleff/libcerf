#include "defs.h" // defines _cerf_cmplx, CMPLX, NaN, Inf
#include <float.h>
#include <math.h>
#include <stdio.h>

typedef struct {
    int failed;
    int total;
} result_t;

// Compute relative error |b-a|/|a|, handling case of NaN and Inf,
static double relerr(double a, double b)
{
    if (isnan(a) || isnan(b) || isinf(a) || isinf(b)) {
        if ((isnan(a) && !isnan(b)) || (!isnan(a) && isnan(b)) || (isinf(a) && !isinf(b))
            || (!isinf(a) && isinf(b)) || (isinf(a) && isinf(b) && a * b < 0))
            return Inf; // "infinite" error
        return 0; // matching infinity/nan results counted as zero error
    }
    if (a == 0)
        return b == 0 ? 0 : Inf;
    else
        return fabs((b - a) / a);
}

// Test whether real numbers 'computed' and 'expected' agree within relative error bound 'limit'
void rtest(result_t* result, double limit, double computed, double expected, const char* name)
{
    ++result->total;
    const double re = relerr(computed, expected);
    if (re > limit) {
        printf("failure in subtest %i: %s\n", result->total, name);
        printf("- fct value %20.15g\n", computed);
        printf("- expected  %20.15g\n", expected);
        printf("=> error %6.2g above limit %6.2g\n", re, limit);
        ++result->failed;
    }
}

// Test whether complex numbers 'computed' and 'expected' agree within relative error bound 'limit'
void ztest(
    result_t* result, double limit, _cerf_cmplx computed, _cerf_cmplx expected, const char* name)
{
    ++result->total;
    const double re_r = relerr(creal(computed), creal(expected));
    const double re_i = relerr(cimag(computed), cimag(expected));
    if (re_r > limit || re_i > limit) {
        printf("failure in subtest %i: %s\n", result->total, name);
        printf("- fct value %20.15g%+20.15g\n", creal(computed), cimag(computed));
        printf("- expected  %20.15g%+20.15g\n", creal(expected), cimag(expected));
        printf("=> error %6.2g or %6.2g above limit %6.2g\n", re_r, re_i, limit);
        ++result->failed;
    }
}

// Wrap rtest; use preprocessor stringification to print the calling 'function_val' verbatim
#define RTEST(result, limit, function_val, expected_val)                                           \
    rtest(&result, limit, function_val, expected_val, #function_val);

// Wrap ztest; use preprocessor stringification to print the calling 'function_val' verbatim
#define ZTEST(result, limit, function_val, expected_val)                                           \
    ztest(&result, limit, function_val, expected_val, #function_val);
