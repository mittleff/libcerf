#!/bin/env python

"""
Computes tau for grid points c_ij and given N, delta.
"""

from mpmath import *
import sys
import derive_w as dw
import enumerate_diameters as ed
import functool as fut
import reltot as rt

if __name__ == '__main__':
    fut.print_provenience()

    N = 20
    delta = 3

    Nb = 16

    nT = int((1.2*Nb)**2) + 1
    D2 = ed.sorted_diameters(nT)

    for ix in range(7*Nb):
        x = ix / Nb
        print(x)

        for iy in range(7*Nb):
            y = iy / Nb
            z = mpc(x, y)
            if abs(z)>7:
                continue
            T = dw.forward(z, 60, False)

            try:
                tau = findroot(lambda tau: rt.rho_of_cNt(z, N, max(0, tau), T) - delta, 0.5,
                               solver='newton', tol=1e-5)
            except Exception as e:
                raise Exception(f'{e} --- z={z}')

            d2 = 0
            for _d2 in reversed(D2):
                if _d2 < (tau*Nb)**2:
                    d2 = _d2
                    break

            print("%10.5g %14.8g %6i" % (y, tau, d2))
        print()
