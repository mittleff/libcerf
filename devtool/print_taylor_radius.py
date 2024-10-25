#!/bin/env python

"""
Compute Taylor coefficients of w(z).
"""

from mpmath import *
import sys
import hp_funcs as hp
import runtool as rt
import functool as fut

mp.dps = 48
mp.pretty = True

def forward(z, N):
    W = []
    W.append(hp.wofz(z.real, z.imag, False))
    W.append(-2*z*W[0] + mpc(0,2)/sqrt(pi))
    for k in range(2,N):
        W.append(-2*(z*W[k-1]+(k-1)*W[k-2]))
        # print("%2i %22.15e" % (k, abs(W[k])))
    R = []
    fac = 1
    for k in range(N):
        R.append(W[k] * fac)
        fac /= (k+1)
    return R

def z2radius(z):
    N = 32
    W = forward(z, N)
    fac = 1
    worst_ratio = 0
    k_worst = 0
    for k in range(1,N):
        ratio = abs(W[k]/W[k-1])
        if ratio > worst_ratio:
            worst_ratio = ratio
            k_worst = k

    radius = 1 / ( 4 * worst_ratio )
    # print("%2i %8g %8g" % (k_worst, worst_ratio, radius))
    return radius

if __name__ == '__main__':
    for y in rt.lingrid(81, 0, 8):
        print(y)
        for x in rt.lingrid(81, 0, 8):
            z = mpc(x, y)
            r = z2radius(z)
            print(x, r)
        print()
