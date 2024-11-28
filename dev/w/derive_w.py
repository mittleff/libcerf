#!/bin/env python

"""
Compute Taylor coefficients of w(z).
"""

from mpmath import *
import sys
sys.path.insert(0, '../shared')
import hp_funcs as hp
import runtool as rt

mp.dps = 48
mp.pretty = True

def dw_integral(z, n):
    return mpc(0,1)**n * sqrt(2)**(n+1) / sqrt(pi) * exp(-z**2/2) * pcfu(n+1/2, mpc(0, -sqrt(2))*z)

def forward(z, N, doublecheck=True):
    """
    Return Taylor coefficients T[k]=w^(k)/k!.
    """
    w2 = mpc(0, -1/sqrt(pi))
    w1 = hp.wofz(z, doublecheck)
    T = [w1]
    for n in range(1, N):
        w1, w2 = -2*(z*w1 + w2)/n, w1
        T.append(w1)
    if N > 0 and doublecheck:
        r = dw_integral(z, N-1)
        if abs(T[-1] - r) > 1e-17*abs(T[-1]):
            print("# %25.16g %25.16g  %25.16g %25.16g  %8g" %
                  (z.real, z.imag, T[-1].real, T[-1].imag, abs(T[-1] - r)/abs(T[-1])))
            raise Exception("w_n: double check failed")
    return T

if __name__ == '__main__':
    if len(sys.argv)==3:
        z = mpc(sys.argv[1], sys.argv[2])
        T = forward(z, 184, False)
        for n in range(len(T)):
            t = T[n]
            r = dw_integral(z, n)
            if abs(t-r)/abs(t)>.5:
                break
            q = 0
            if n>1:
                q = abs(T[n]/T[n-1])
            print("%2i %11.3e %11.3e %8e %8e %8.3f" %
                  (n, t.real, t.imag, abs(t), abs(t-r)/abs(t), q) )
        sys.exit(0)

    N = 24
    if len(sys.argv)==2:
        N = int(sys.argv[1])

    for x in [0] + rt.loggrid(20, 1e-4, 1.):
        print(x)
        for y in [0] + rt.loggrid(20, 1e-4, 1.):
            z = mpc(x, y)
            T = forward(z, N)
            t=T[-1]
            r = dw_integral(z, N-1)
            print("%13.5e %11.3e %11.3e %11.3e %11.3e %11.3e %8g" %
                  (y, t.real, t.imag, abs(t), r.real, r.imag, abs(t-r)/abs(t)))
        print()
