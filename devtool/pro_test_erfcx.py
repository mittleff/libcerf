#!/bin/env python

"""
Compute Chebyshev coefficients for erfcx.
Used to generate code sections in erfcx.c.
"""

from mpmath import *
import random, sys

mp.dps = 48
mp.pretty = True

tuning = True
final = False # Extra checks, to be turned on in final production run

def dawson_kernel(t):
    return exp(t**2)

def hp_erfcx(x):
    result = exp(x**2)*erfc(x)
    if False: # TODO restore, currently copied&pasted from im_wofx: final:
        # Check mpmath-computed reference value against mpmath-based brute-force integration
        r2 = 2/sqrt(pi) * exp(-x**2) * quad(dawson_kernel, [0, x])
        if abs(result-r2)/result > 1e-17:
            raise Exception(f"mpmath inaccurate")
    return result

def tol(x):
    if x < -16:
        return 2e-13
    if x < -8:
        return 8e-14
    if x < -4:
        return 1.6e-14
    if x < -2:
        return 3e-15
    if x < -1:
        return 9e-16
    return 4.5e-16

if __name__ == '__main__':

    X = []

    Ni = 957
    x_fr = .06
    x_to = 200
    step = log10(x_to/x_fr)/(Ni-1)
    for i in range(Ni):
        x = x_fr * 10**(i*step)
        if x < 26.6:
            X.append(-x)
        X.append(x)

    Ni = 357
    x_fr = 200
    x_to = 2e8
    step = log10(x_to/x_fr)/(Ni-1)
    for i in range(Ni):
        x = x_fr * 10**(i*step)
        X.append(x)

    print("//--- The following code is generated by " + sys.argv[0])
    print("// clang-format off")

    for x in sorted(X):
        xs = "%+23.15e" % x
        print("    RTEST(result, %g, erfcx(%s), %+23.16e);" % (tol(x), xs, hp_erfcx(mpf(xs))))

    print("// clang-format on")
    print("//--- End of autogenerated code")
