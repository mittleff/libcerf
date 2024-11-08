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

def dw_integral(z, n):
    return mpc(0,1)**n * sqrt(2)**(n+1) / sqrt(pi) * exp(-z**2/2) * pcfu(n+1/2, mpc(0, -sqrt(2))*z)

def forward(z, N):
    """
    Return Taylor coefficients T[k]=w^(k)/k!.
    """
    w2 = mpc(0, -1/sqrt(pi))
    w1 = hp.wofz(z, False)
    T = [w1]
    for n in range(1, N):
        w1, w2 = -2*(z*w1 + w2)/n, w1
        T.append(w1)
    if abs(T[-1] - dw_integral(z, N-1)) > 1e-17*abs(T[-1]):
        raise Exception("w_n: double check failed")
    return T

if __name__ == '__main__':
    if len(sys.argv)==3:
        z = mpc(sys.argv[1], sys.argv[2])
        T = forward(z, 60)
        for n in range(len(T)):
            t = T[n]
            r = dw_integral(z, n)
            q = 0
            p = 0
            if n>0:
                p = abs(T[n])/abs(T[n-1])
                q = abs(R[n])/abs(R[n-1])
            print("%11.3e %11.3e %11.3e %8f  %11.3e %11.3e %11.3e %8f  %11.3e %11.3e" %
                  (t.real, t.imag, abs(t), p, r.real, r.imag, abs(r), q, B[n].real, B[n].imag), )
        sys.exit(0)

    N = 24
    for x in [0] + rt.loggrid(30, 1e-18, .08) + rt.loggrid(50, .1, 7.):
        print(x)
        for y in [0] + rt.loggrid(30, 1e-18, .08) + rt.loggrid(50, .1, 7.):
            z = mpc(x, y)
            if z==0:
                continue
            T = forward(z, N)
            t=T[-1]
            r = dw_integral(z, N-1)
            print("%13.5e %11.3e %11.3e %11.3e %11.3e %8g" %
                  (y, t.real, t.imag, r.real, r.imag, abs(t-r)/abs(t)))
        print()
