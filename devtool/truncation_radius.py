#!/bin/env python

"""
For given order N and central z, print RTE(tau).
"""

from mpmath import *
import sys
import hp_funcs as hp
import functool as fut
import derive_w as dw
import runtool as rt

mp.dps = 48
mp.pretty = True

def taylor_remainder(z, t, N, M, T):
    r = 0
    for n in range(N, M):
        r += abs(T[n]) * t**n
    return r

def rte(z, t, N):
    amin = abs(hp.wofz(z+1/sqrt(2)*mpc(1,1)))

    T = dw.forward(z, 100, False)
    amax = taylor_remainder(z, t, N, 90, T)
    amax2 = taylor_remainder(z, t, N, 100, T)
    if abs((amax2-amax)/amax) > 1e-17:
        raise Exception("Truncation error in computation of TE of w(z)")

    xi = amax / amin

    return xi

if __name__ == '__main__':
    if len(sys.argv)==4:
        N = int(sys.argv[1])
        z = mpc(sys.argv[2], sys.argv[3])
        for t in rt.loggrid(40, .01, 2):
            xi = rte(z, t, N)
            print("%12.5g %12.5g" % (t, xi))
