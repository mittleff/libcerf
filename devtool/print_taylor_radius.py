#!/bin/env python

"""
Compute Taylor coefficients of w(z).
"""

from mpmath import *
import sys
import hp_funcs as hp
import functool as fut
import derive_w as dw

mp.dps = 48
mp.pretty = True

def N2rho(N):
    f = lambda r: 2**-53*(1-2*r) - r**N
    rho = findroot(f, [.01, .49], solver='bisect')
    return rho

def z2radius(z, N, rho):
    """
    For a Taylor expansion around z, determine radius R such that T[k]*R/T[k-1] < rho.
    """

    T = dw.forward(z, 3*N//2)
    fac = 1
    phi = 0
    k_worst = 0
    for k in range(1, len(T)):
        ratio = abs(T[k]/T[k-1])
        if ratio > phi:
            phi = ratio
            k_worst = k
    if k_worst >= N:
        for k in range(len(T)):
            t = T[k]
            q = 0
            if k>0:
                q = abs(t)/abs(T[k-1])
            print("#>> %3i %10e %10e %8f" % (k, t.real, t.imag, q))
        raise Exception("bound phi_N is insufficient; f^(n) don't decay fast enough")
    R = rho / phi

    sum = abs(T[0])
    Nused = len(T)
    for k in range(1,len(T)):
        term = abs(T[k]) * R**k
        if term < 2**-53 * sum:
            Nused = k-1
            break
        sum -= term
    if Nused > N:
        raise Exception("determination of R failed")

    return R

if __name__ == '__main__':
    N = 24
    if len(sys.argv)==3:
        z = mpc(sys.argv[1], sys.argv[2])
        R = z2radius(z, N)
        e = int(100*R/sqrt(2))/100
        print("N R %8g d^2 %8g edge %5.2f x range %5.2f %5.2f y range %5.2f %5.2f" %
              (R, 4*R**2, e, z.real-e, z.real+e, z.imag-e, z.imag+e))
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
            R = z2radius(z, N)
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
