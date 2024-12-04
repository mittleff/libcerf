#!/bin/env python

"""
Computes tau for grid points c_ij and given N, delta.
"""

from mpmath import *
import derive_w as dw
import sys
sys.path.insert(0, '../shared')
import hp_funcs as hp
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

def wmin(ix, iy, Nb, d2):
    # print(f"WMIN ix={ix} iy={iy} d2={d2}")
    wmin = inf
    for dix in range(ix%2, int(sqrt(d2)), 2):
        diy = iy%2 + int((sqrt(d2-dix**2)-iy%2)/2)
        wmin = min(wmin, abs(hp.wofz(mpc((ix+dix)/Nb, (iy+diy)/Nb))))
        #print(f"... dix={dix} diy={diy} wmin={'%12g' % wmin}")
    #print(f"RETURN wmin={'%12g' % wmin}")
    return wmin

if __name__ == '__main__':
    if len(sys.argv)!=3:
        print(f"Usage: {sys.argv[0]} N delta")
        sys.exit(1)
    N = int(sys.argv[1])
    delta = float(sys.argv[2])

    fut.print_provenience()

    d2max = 400
    S = [[sorted_diameters(d2max, sx, sy) for sy in [0, 1]] for sx in [0, 1]]

    Nb = 16

    nT = int((1.2*Nb)**2) + 1

    for ix in range(7*Nb+1):
        x = ix / Nb
        print(x)

        for iy in range(7*Nb+1):
            y = iy / Nb
            z = mpc(x, y)
            if abs(z)>7:
                continue
            Sxy = S[ix%2][iy%2]
            T = dw.forward(z, 60, False)
            assert(rt.rho_of_cNt(z, N, 0, T) <= delta)

            # Binary search for maximum itau such that relerr <= delta.
            nGenS = int(math.log(len(Sxy), 2)) # number of binary generations
            itau = 0
            for n in reversed(range(nGenS)):
                inew = itau + 2**n
                if inew >= len(Sxy):
                    continue
                d2 = Sxy[inew]
                tau = sqrt(d2) / Nb
                if rt.abserr(z, N, tau, T) / wmin(ix, iy, Nb, d2) <= delta:
                    itau = inew

            if itau >= len(Sxy)-1:
                raise Exception("increase d2max!")
            print("%10.5g %3i" % (y, Sxy[itau]))
        print()
