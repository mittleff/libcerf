#!/bin/env python

"""
Compute Chebyshev coefficients for erfcx.
Used to generate code sections in erfcx.c.
"""

from mpmath import *
import functool as fut
import random, sys

mp.dps = 48
mp.pretty = True

tuning = True
final = False # Extra checks, to be turned on in final production run

def dawson_kernel(t):
    return exp(t**2)

def highprecision_erfcx(x):
    result = exp(x**2)*erfc(x)
    if False: # TODO restore, currently copied&pasted from im_wofx: final:
        # Check mpmath-computed reference value against mpmath-based brute-force integration
        r2 = 2/sqrt(pi) * exp(-x**2) * quad(dawson_kernel, [0, x])
        if abs(result-r2)/result > 1e-17:
            raise Exception(f"mpmath inaccurate")
    return result

if __name__ == '__main__':

    Nout = 9

    n2e = 6
    R = []

    nOctaves = 7
    iOctave0 = -3
    for ir in range(nOctaves):
        for js in range(2**n2e):
            a = (2**n2e+js) * 2**(ir-n2e+iOctave0)
            b = a + 2.**(ir-n2e+iOctave0)
            # if b>=.5 and a<=12.:
            R.append((a,b, ir, js))
    C = fut.chebcoef(R, Nout, highprecision_erfcx)

    # fut.print_cheby_coeffs(C)
    # fut.print_clenshaw_code(R, C, Nout)
    fut.print_powerseries_code(R, C, Nout)
