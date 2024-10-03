#!/bin/env python

"""
Computes accuracy of external implementation of function w_of_z
by comparison with an internal high-precision function.

The comparison is done for z=r2z(r,s) with r,s from some grids,
and for additional points obtained by bisection around points r
where the employed algorithm or number of terms has changed.

The grids are hard-coded here, and are changed ad hoc during development.

Output is either a table with one line for each r,
or a report on the worst case (the r with largest inaccuracy of f(x)).
"""

from mpmath import *
import sys
import hp_funcs as hp

mp.dps = 48
mp.pretty = True

if __name__ == '__main__':
    c = mpc(5.5, .5)
    z = mpc(5.9, .9)
    wz = hp.wofz(z.real, z.imag)
    d = z - c
    fac = 1
    sum = 0
    for k in range(0,20):
        if k==0:
            w0 = hp.wofz(c.real, c.imag)
            wk = w0
        elif k==1:
            w1 = -2*c*w0 + mpc(0,2)/sqrt(pi)
            wk = w1
        else:
            wk = -2*(c*w1+(k-1)*w0)
            w0 = w1
            w1 = wk
            fac /= k
        sum += fac*wk*d**k
        print("%2i %8e %22.15e %22.15e  %8e" %
              (k, fac*abs(wk), sum.real, sum.imag, abs(sum-wz)))
