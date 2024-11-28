#!/bin/env python

"""
For given order N and expansion center c, print relative error rho as function of tau.
This script was used to generate figure cgt/fig/reltot.eps.
"""

from mpmath import *
import sys
sys.path.insert(0, '../shared')
import hp_funcs as hp
import functool as fut
import derive_w as dw
import runtool as rt

mp.dps = 48
mp.pretty = True

def taylor_remainder(t, N, M, T):
    r = 0
    for n in range(N, M):
        r += abs(T[n]) * t**n
    return r

def rho_of_cNt(z, N, tau, T):
    """
    Returns total relative error in units of epsilon.
    """
    eps = 2**-53

    # Minimum function value (denominator):
    amin = abs(hp.wofz(z+tau/sqrt(2)*mpc(1,1)))

    # Truncation error:
    dt = taylor_remainder(tau, N, 50, T)
    ddt = taylor_remainder(tau, 50, 60, T)
    if abs(ddt) > 1e-13 * dt:
        if dt>1e3:
            dt = inf
        if dt<1e-2 and ddt<1e-4:
            dt = 0
        else:
            raise Exception("Truncation error in computation of TE of w(z): N=%i tau=%g, dt=%24.16g, ddt=%24.16g" % (N, tau, dt, ddt))
    te = dt / eps

    # Rounding error:
    la = sqrt(5)
    if z.imag==0 or z.real==0:
        la = 1
    re = 0
    for n in range(N):
        vn = hypot(float(T[n].real)-T[n].real, float(T[n].imag)-T[n].imag)/abs(T[n])/eps
        cn = 2
        if n==0 or n==N-1:
            cn = 1
        re += abs(T[n]) * tau**n * (vn + cn + n*(1+la))

    # Total relative error:
    return (te + re) /amin

if __name__ == '__main__':
    if len(sys.argv)==3:
        z = mpc(sys.argv[1], sys.argv[2])
        la = sqrt(5)
        if z.imag==0 or z.real==0:
            la = 1
        tmax = 3
        Nt = 780
        T = dw.forward(z, 60, False)
        for N in (16,20,24,28,32):
            print(N)
            for tau in rt.lingrid(Nt, tmax/Nt, tmax):
                rho = rho_of_cNt(z, N, tau, T)
                print("%12.5g %12.5g" % (tau, rho))
            print()
