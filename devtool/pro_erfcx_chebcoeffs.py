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

final = False # Extra checks, to be turned on in final production run

def hp_erfcx(x, doublecheck=False):
    result = exp(x**2)*erfc(x)
    if doublecheck:
        # Check mpmath-computed reference value against mpmath-based brute-force integration
        if x<6:
            r2 = exp(x**2) - 2/sqrt(pi) * quad(lambda t : exp(x**2-t**2), [0, x])
        else:
            r2 = 2/sqrt(pi) * quad(lambda t : exp(x**2-t**2), [x, mpf('inf')])
        if abs(result-r2)/result > 1e-17:
            raise Exception(f"mpmath inaccurate")
    return result

if __name__ == '__main__':

    Nout = 9

    R = fut.octavicRanges(.125, 16, 6)
    C = fut.chebcoef(R, Nout, hp_erfcx, final)

    # fut.print_cheby_coeffs(C)
    # fut.print_clenshaw_code(R, C, Nout)
    fut.print_powerseries_code(R, C, Nout)
