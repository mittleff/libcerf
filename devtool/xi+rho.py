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

mp.dps = 80
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

                amax = taylor_remainder(z, tau, N, 90, T)
                amax2 = taylor_remainder(z, tau, N, 100, T)
                if abs((amax2-amax)/amax) > 1e-13:
                    if amax>1e3 and amax2>1e3:
                        amax = inf
                    if amax<1e-2 and amax2<1e-2:
                        amax = 0
                    else:
                        raise Exception("Truncation error in computation of TE of w(z): N=%i tau=%g, amax=%24.16g, amax2=%24.16g" % (N, tau, amax, amax2))
                xi = amax / amin / eps

                re = 0   # using actual rounding errors vn
                rewc = 0 # wc = worst case
                for n in range(N):
                    vn = hypot(float(T[n].real)-T[n].real, float(T[n].imag)-T[n].imag)/abs(T[n])/eps
                    vnwc = 1
                    if z==0 and n==0:
                        vnwc = 0
                    cn = 2
                    if n==0 or n==N-1:
                        cn = 1
                    re += abs(T[n]) * tau**n * (vn + cn + n*(1+la))
                    rewc += abs(T[n]) * tau**n * (vnwc + cn + n*(1+la))
                rho = re/amin
                rhowc = rewc/amin
                print("%12.5g %12.5g %12.5g %12.5g %12.5g %12.5g" % (tau, xi, rho, rhowc, xi+rho, xi+rhowc))
            print()
