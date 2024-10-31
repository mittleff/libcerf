#!/bin/env python

"""
Compute Taylor coefficients of w(z).
"""

from mpmath import *
import sys, re
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

def check_zr(z, R):
    """
    """
    N = 36
    T = forward(z, N)
    fac = 1
    worst_ratio = 0
    for k in range(1,N):
        ratio = abs(T[k]*R/T[k-1])
        if ratio > worst_ratio:
            worst_ratio = ratio

    sum = abs(T[0])
    Nmax = 99
    for k in range(1,len(T)):
        term = abs(T[k]) * R**k
        if term < 2**-54 * sum:
            Nmax = k-1
            break
        sum -= term

    return worst_ratio, Nmax

def check(ix, iy, ir):
    z = mpc(ix * 0.125, iy * 0.125)
    R = ir * 0.125 * sqrt(2)
    worst_ratio, Nmax = check_zr(z, R)
    warn = ""
    limit = .25
    if worst_ratio > limit:
        warn = "CANCELLATION"
    elif Nmax>22:
        warn = "N"
    print("x %2i y %2i N %2i w' %5.3f %s" % (ix, iy, Nmax, worst_ratio/limit, warn))

if __name__ == '__main__':
    with open('devtool/tiling.ps', 'r') as file:
        for line in file:
            m = re.match(r'^\s*(\d+)\s*(\d+)\s*(\d+) Q\s*$', line)
            if m:
                ix = int(m[1])
                iy = int(m[2])
                ir = int(m[3])
                check(ix, iy, ir)
