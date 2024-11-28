#!/bin/env python

"""
Print expansion coefficients of erfcx for the non-Chebyshev cases.
"""

from mpmath import *

mp.dps = 64
mp.pretty = True

eps = 2**(-53)
I = mpc(0,1)

if __name__ == '__main__':
    # Taylor: decaying too slowly, not investigated further.
    for n in reversed(range(0,15)):
        c = pow(-1,n) / gamma(n/2.+1)
        print("%+23.16e ) * x" % c)
