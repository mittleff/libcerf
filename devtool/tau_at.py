#!/bin/env python

"""
Computes tau for given c, N, delta.
"""

from mpmath import *
import sys
import reltot as rt

if __name__ == '__main__':
    if len(sys.argv)!=3:
        raise Exception(f'Usage: {sys.argv[0]} zre zim')

    z = mpc(sys.argv[1], sys.argv[2])

    N = 20
    delta = 3

    tau = findroot(lambda tau: rt.rho_of_cNt(z, N, tau) - delta, 0.1, solver='newton')
    print("%12.5g" % (tau))
