#!/bin/env python

"""
Compute Chebyshev coefficients for im_w_of_x.
Used to generate code sections in im_w_of_z.c.
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

def highprecision_imwx(x):
    fz = exp(-x**2)*erfc(mpc(0, -x))
    result = fz.imag
    if final:
        # Check mpmath-computed reference value against mpmath-based brute-force integration
        r2 = 2/sqrt(pi) * exp(-x**2) * quad(dawson_kernel, [0, x])
        if abs(result-r2)/result > 1e-17:
            raise Exception(f"mpmath inaccurate")
    return fz.imag

if __name__ == '__main__':

    Nout = 8

    n2e = 6
    R = []
    for ir in range(5):
        for js in range(2**n2e):
            a = (2**n2e+js) * 2**(ir-n2e-1)
            b = a + 2.**(ir-n2e-1)
            if b>=.5 and a<=12.:
                R.append((a,b, ir, js))
    C = fut.chebcoef(R, Nout, highprecision_imwx)

    # fut.print_cheby_coeffs(C)
    # fut.print_clenshaw_code(R, C, Nout)
    fut.print_powerseries_code(R, C, Nout)
