#!/bin/env python

"""
Compute Taylor coefficients of w(z).
"""

from mpmath import *
import sys
import hp_funcs as hp

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
    w2, w1 = 1/gamma(N+dN)**.67, 0
    recurse = lambda n, w1, w2: (-z*w1 - n*w2/2, w1)
    for n in reversed(range(N, N+dN)):
        w1, w2 = recurse(n, w1, w2)
    B = []
    for n in reversed(range(N)):
        w1, w2 = recurse(n, w1, w2)
        B.append(w1)
    fac = hp.wofz(z.real, z.imag, True) / B[-1]
    R = []
    for b in reversed(B):
        R.append(fac * b)
    return R

if __name__ == '__main__':
    if len(sys.argv)==3:
        z = mpc(sys.argv[1], sys.argv[2])
        T = forward(z, 60)
        R = backward(z, 60, 60)
        for n in range(len(T)):
            t = T[n]
            r = R[n]
            q = 0
            p = 0
            if n>0:
                p = abs(T[n])/abs(T[n-1])
                q = abs(R[n])/abs(R[n-1])
            print("%11.3e %11.3e %11.3e %8f  %11.3e %11.3e %11.3e %8f" %
                  (t.real, t.imag, abs(t), p, r.real, r.imag, abs(r), q))
        sys.exit(0)
