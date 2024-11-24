#!/bin/env python

"""
For given order N and central z, print functions xi (RTE) and rho (RRE).
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

    return xi

if __name__ == '__main__':
    eps = 2**-53
    if len(sys.argv)==3:
        z = mpc(sys.argv[1], sys.argv[2])
        la = sqrt(5)
        if z.imag==0 or z.real==0:
            la = 1
        tmax = 3
        Nt = 780
        for N in (16,20,24,28,32):
            print(N)
            for tau in rt.lingrid(Nt, tmax/Nt, tmax):

                T = dw.forward(z, 120, False)

                amin = abs(hp.wofz(z+tau/sqrt(2)*mpc(1,1)))

                amax = taylor_remainder(z, tau, N, 50, T)
                amax2 = taylor_remainder(z, tau, N, 60, T)
                if abs((amax2-amax)/amax) > 1e-17:
                    raise Exception("Truncation error in computation of TE of w(z)")
                xi = amax / amin / eps

                re = 0
                for n in range(N):
                    re += abs(T[n]) * tau**n * (n + 1 + n*la)
                rho = re/amin

                print("%12.5g %12.5g %12.5g" % (tau, xi, rho))
            print()
