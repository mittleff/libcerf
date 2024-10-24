#!/bin/env python

"""
Data for plot
"""

from mpmath import *
import sys
import hp_funcs as hp

mp.dps = 48
mp.pretty = True

if __name__ == '__main__':
    if len(sys.argv)<4:
        print(f"Usage: {sys.argv[0]} M Re(z) Im(z)")
        sys.exit(0)
    M = int(sys.argv[1])
    z = mpc(sys.argv[2], sys.argv[3])

    rho = z
    for k in reversed(range(1,M)):
        rho = z - k/2/rho
        print("%2i %22.15e %22.15e" % (k, rho.real, rho.imag))

    w = mpc(0,1)/sqrt(pi)/rho
    wref = hp.wofz(z.real, z.imag, True)
    if not hp.cagree(w, wref):
        raise Exception(f"inaccurate")
