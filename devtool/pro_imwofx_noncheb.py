#!/bin/env python

"""
Print expansion coefficients of Im w(x) for the non-Chebyshev cases.
"""

from mpmath import *

mp.dps = 64
mp.pretty = True

eps = 2**(-53)
I = mpc(0,1)

fodd = [1] # factorial (2n-1)!!
for n in range(1,20):
    fodd.append((2*n-1) * fodd[-1])

if __name__ == '__main__':
    for n in reversed(range(0,11)):
        c = fodd[n] / 2**n / sqrt(pi)
        print(" + %22.16e ) * (r*r)" % c)
