#!/bin/env python

"""
Compute Taylor coefficients of w(z).
"""

from mpmath import *
import sys
import hp_funcs as hp
import runtool as rt

mp.dps = 48
mp.pretty = True

def forward(z):
    W = []
    W.append(hp.wofz(z.real, z.imag, True))
    W.append(-2*z*W[0] + mpc(0,2)/sqrt(pi))
    fac = 1
    for k in range(2,20):
        fac /= k
        W.append(fac * (-2*(z*W[k-1]+(k-1)*W[k-2])))
        # print("%2i %22.15e" % (k, abs(W[k])))
    return W

def z2radius(z):
    W = forward(z)
    fac = 1
    worst_ratio = 0
    k_worst = 0
    for k in range(1,20):
        ratio = abs(W[k]/W[k-1])
        if ratio > worst_ratio:
            worst_ratio = ratio
            k_worst = k

    radius = 1 / ( 4 * worst_ratio )
    # print("%2i %8g %8g" % (k_worst, worst_ratio, radius))
    return radius

if __name__ == '__main__':
    # if len(sys.argv)<3:
    #     print(f"Usage: {sys.argv[0]} Re(z) Im(z)")
    #     sys.exit(0)

    for y in rt.loggrid(201, .1, 10):
        print(y)
        for x in rt.loggrid(201, .1, 10):
            z = mpc(x, y)
            r = z2radius(z)
            print(x, r)
        print()
