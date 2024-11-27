#!/bin/env python

"""
Computes tau for grid points c_ij and given N, delta.
"""

from mpmath import *
import sys
import reltot as rt

if __name__ == '__main__':

    N = 20
    delta = 3

    Nb = 16
    for ix in range(7*Nb):
        x = ix / Nb
        print(x)

        for iy in range(7*Nb):
            y = iy / Nb
            z = mpc(x, y)
            if abs(z)>7:
                continue

            tau = findroot(lambda tau: rt.rho_of_cNt(z, N, tau) - delta, 0.1,
                           solver='newton', tol=1e-5)
            print("%12.5g %12.5g" % (y, tau))
        print()
