#!/bin/env python

"""
Computes tau for grid points c_ij and given N, delta.
"""

from mpmath import *
import derive_w as dw
import enumerate_diameters as ed
import sys
sys.path.insert(0, '../shared')
import functool as fut
import reltot as rt
import bisect, math

def sorted_diameters(N, sx, sy):
    D2 = set()
    for j in range(1+sy, N, 2):
        for i in range(1+sx, N, 2):
            d2 = j**2 + i**2
            if d2 > N:
                continue
            D2.add(d2)
    return sorted(D2)


if __name__ == '__main__':
    fut.print_provenience()

    d2max = 400
    S = [[sorted_diameters(d2max, sx, sy) for sy in [0, 1]] for sx in [0, 1]]

    N = 20
    delta = 3

    Nb = 16

    nT = int((1.2*Nb)**2) + 1
    D2 = ed.sorted_diameters(nT)

    for ix in range(7*Nb+1):
        x = ix / Nb
        print(x)

        for iy in range(7*Nb+1):
            y = iy / Nb
            z = mpc(x, y)
            if abs(z)>7:
                continue

            T = dw.forward(z, 60, False)
            assert(rt.rho_of_cNt(z, N, 0, T) <= delta)

            Sxy = S[ix%2][iy%2]
            nGenS = int(math.log(len(Sxy), 2)) # number of binary generations
            itau = 0
            for n in reversed(range(nGenS)):
                inew = itau + 2**n
                if inew >= len(Sxy):
                    continue
                tau = sqrt(Sxy[inew]) / Nb
                if rt.rho_of_cNt(z, N, tau, T) <= delta:
                    itau = inew

            print("%10.5g %8.5g %3i %3i" % (y, tau, Sxy[itau], itau))
        print()
