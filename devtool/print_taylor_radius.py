#!/bin/env python

"""
Compute Taylor coefficients of w(z).
"""

from mpmath import *
import sys
import hp_funcs as hp
import runtool as rt
import functool as fut
import enumerate_polyominoes as ep

mp.dps = 48
mp.pretty = True

Nmax = 24
f = lambda r: 2**-53*(1-2*r) - r**Nmax
rho = findroot(f, [.01, .49], solver='bisect')

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

def errN(T, R):
    """
    Return maximum relative error of T[N-1] term in Taylor sum.
    """
    sum = abs(T[0])
    for k in range(1,Nmax):
        term = abs(T[k]) * R**k
        sum -= term
    tN = abs(T[Nmax]) * R**Nmax
    return tN/sum

def redetermine_R(T, R):
    """
    Called if Nused > N. Recomputes R such that Nused = N.
    """
    f = lambda r: errN(T, r) - 2**-53
    r = mp.findroot(f, [R], solver='newton', verbose=False)
    return r, Nmax

def z2radius(z):
    """
    For a Taylor expansion around z, determine radius R such that T[k]*R/T[k-1] < q,
    thereby ensuring that sum T[k] >= (1-2q)/(1-q) T[0]
    """
    Ntest = Nmax+2
    T = forward(z, Ntest)
    fac = 1
    worst_ratio = 0
    k_worst = 0
    for k in range(1,Ntest):
        ratio = abs(T[k]/T[k-1])
        if ratio > worst_ratio:
            worst_ratio = ratio
            k_worst = k
    R = rho / worst_ratio

    sum = abs(T[0])
    Nused = Ntest
    for k in range(1,len(T)):
        term = abs(T[k]) * R**k
        if term < 2**-53 * sum:
            Nused = k-1
            break
        sum -= term

    if Nused > Nmax:
        # print("original:", R, Nmax)
        R, Nused = redetermine_R(T, R)
    return R, Nused

if __name__ == '__main__':
    if len(sys.argv)==3:
        z = mpc(sys.argv[1], sys.argv[2])
        R, Nmax = z2radius(z)
        e = int(100*R/sqrt(2))/100
        print("N %2i R %8g d^2 %8g edge %5.2f x range %5.2f %5.2f y range %5.2f %5.2f" %
              (Nmax, R, 4*R**2, e, z.real-e, z.real+e, z.imag-e, z.imag+e))
        sys.exit(0)

    P = ep.sorted_polyominoes(100)
    Diams2 = []
    Rows = []
    for i, j, d2, area, rows in P:
        Diams2.append(d2)
        Rows.append(rows)

    fut.print_provenience()
    for iy in range(128):
        y = iy/16
        # print(iy)
        for ix in range(128):
            x = ix/16
            z = mpc(x, y)
            R, Nmax = z2radius(z)
            ir = 997
            for ii in range(len(Diams2)):
                if (2*8*R)**2 < Diams2[ii]:
                    ir = ii - 1
                    break
            assert(ir < 997)
            while (iy%2 != len(Rows[ir])%2 or ix%2 != Rows[ir][0]%2) and ir >= 0:
                ir -= 1
            print("(%3i, %3i, %3i, %3i, %8.4f), # %s" % (iy, ix, Nmax, ir, R, str(Rows[ir])))
        # print()
