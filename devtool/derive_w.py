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

def forward(z, N):
    """
    Return Taylor coefficients T[k]=w^(k)/k!.
    """
    W = []
    W.append(hp.wofz(z.real, z.imag, False))
    W.append(-2*z*W[0] + mpc(0,2)/sqrt(pi))
    for k in range(2,N):
        W.append(-2*(z*W[k-1]+(k-1)*W[k-2]))
        # print("%2i %22.15e" % (k, abs(W[k])))
    T = []
    fac = 1
    for k in range(N):
        T.append(W[k] * fac)
        fac /= (k+1)
    return T

def backward(z, N, dN):
    """
    Return Taylor coefficients T[k]=w^(k)/k!,
    computed through backward recursion.
    """
    recurse = lambda n, w1, w2: (-z*w1 - (n+2)*w2/2, w1)
    w2, w1 = mpc(2,1)/gamma(N+dN)**.67, 0
    for n in reversed(range(N, N+dN)):
        w1, w2 = recurse(n, w1, w2)

    # check accuracy
    v2, v1 = mpc(-1,-2)/gamma(N+dN+7)**.67, 0
    for n in reversed(range(N, N+dN+7)):
        v1, v2 = recurse(n, v1, v2)
    if abs(w1/w2 - v1/v2) > 1e-16 * (abs(w1/w2) + abs(v1/v2)):
        pass #raise Exception("dN too small?")

    B = []
    for n in reversed(range(N)):
        w1, w2 = recurse(n, w1, w2)
        B.insert(0, w1)
    w1, w2 = recurse(-1, w1, w2)
    fac = mpc(0, -1/sqrt(pi)) / w1

    R = [fac * b for b in B]

    return R, B

def dw_integral(z, n):
    return 2/sqrt(pi) * 2**n/gamma(n+1) * quad(lambda u: u**n * exp(-u**2-mpc(0,2)*u*z), [0, +inf])

if __name__ == '__main__':
    if len(sys.argv)==3:
        z = mpc(sys.argv[1], sys.argv[2])
        T = forward(z, 60)
        R, B = backward(z, 60, 40)
        for n in range(len(T)):
            t = T[n]
            r = R[n]
            q = 0
            p = 0
            if n>0:
                p = abs(T[n])/abs(T[n-1])
                q = abs(R[n])/abs(R[n-1])
            print("%11.3e %11.3e %11.3e %8f  %11.3e %11.3e %11.3e %8f  %11.3e %11.3e" %
                  (t.real, t.imag, abs(t), p, r.real, r.imag, abs(r), q, B[n].real, B[n].imag), )
        sys.exit(0)

    for x in [0] + rt.loggrid(30, 1e-18, .08) + rt.loggrid(50, .1, 7.):
        print(x)
        for y in [0] + rt.loggrid(30, 1e-18, .08) + rt.loggrid(50, .1, 7.):
            z = mpc(x, y)
            if z==0:
                continue
            T = forward(z, 25)
            t=T[24]
            r = dw_integral(z, 24)
            print("%13.5e %11.3e %11.3e %11.3e %11.3e %8g" %
                  (y, t.real, t.imag, r.real, r.imag, abs(t-r)/abs(t)))
        print()
