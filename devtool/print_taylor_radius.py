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

def err23(T, R):
    """
    Return maximum relative error of T[23] term in Taylor sum.
    """
    sum = abs(T[0])
    for k in range(1,22):
        term = abs(T[k]) * R**k
        sum -= term
    t23 = abs(T[23]) * R**23
    # print(f"R {R} t {t23} s {sum} a {t23/sum - 2**-54}")
    return t23/sum

def redetermine_R(T, R):
    """
    Called if Nmax > 23. Recomputes R such that Nmax = 23.
    """
    f = lambda r: err23(T, r) - 2**-54
    # print(f"1.0: {R} -> {f(R)}")
    # print(f"0.8: {0.8*R} -> {f(0.8*R)}")
    r = mp.findroot(f, [R], solver='newton', verbose=False)
    return r, 23

def z2radius(z):
    """
    For a Taylor expansion around z, determine radius R such that T[k]*R/T[k-1] < 1/4,
    thereby ensuring that sum T[k] >= 2/3 T[0]
    """
    N = 32
    T = forward(z, N)
    fac = 1
    worst_ratio = 0
    k_worst = 0
    for k in range(1,N):
        ratio = abs(T[k]/T[k-1])
        if ratio > worst_ratio:
            worst_ratio = ratio
            k_worst = k
    R = 1 / ( 4 * worst_ratio )

    sum = abs(T[0])
    Nmax = 99
    for k in range(1,len(T)):
        term = abs(T[k]) * R**k
        if term < 2**-54 * sum:
            Nmax = k-1
            break
        sum -= term

    if Nmax > 23:
        print("original:", R, Nmax)
        R, Nmax = redetermine_R(T, R)
    return R, Nmax

def zr2nmax(z, R):
    return 99

if __name__ == '__main__':
    if len(sys.argv)==3:
        z = mpc(sys.argv[1], sys.argv[2])
        R, Nmax = z2radius(z)
        print(R, Nmax)
        sys.exit(0)

    for y in rt.lingrid(81, 0, 8):
        print(y)
        for x in rt.lingrid(81, 0, 8):
            z = mpc(x, y)
            R, Nmax = z2radius(z)
            print(x, R, Nmax)
        print()
