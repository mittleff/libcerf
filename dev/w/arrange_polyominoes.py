#!/bin/env python

"""
Generates one possible polyomino coverage.
"""

import sys
sys.path.insert(0, '../shared')
from mpmath import *
from math import sqrt
import hp_funcs as hp
import functool as fut
import enumerate_polyominoes as ep
import datetime

mp.dps = 48
mp.pretty = True

def sorted_diameters(N, sx, sy):
    """
    Returns sorted list of possible values of squared circumcircle diameter
    in units of a, for given parities of center.
    """
    D2 = set()
    for j in range(1+sy, N, 2):
        for i in range(1+sx, N, 2):
            d2 = j**2 + i**2
            if d2 > N:
                continue
            D2.add(d2)
    return sorted(D2)

def polyomino_pattern(d2, sx, sy):
    rows = []
    nj2 = int(sqrt(d2)) + 1
    for k in range(2-sx, nj2, 2):
        dx = sqrt(d2 - k**2)/2 - (sy)/2
        row = 2 * int(dx) + (sy)
        if row > 0:
            if k > sx:
                rows.insert(0, row)
            rows.append(row)
    return rows

def read_d2_file(fname):
    """
    Reads file that provides d2(x,y).
    Returns Nb=1/b and D2 (squared diameters for b-lattice points).
    """
    with open(fname, 'r') as f:
        t = ''
        for line in f:
            if not line[0] =='#':
                t += line
    Nb = None
    D2 = []
    for block in t.split('\n\n'):
        if block=='':
            break
        a = block.split('\n')
        x = float(a[0])
        if x == 0:
            ix = 0
        else:
            ix += 1
            if round(Nb*x) != ix:
                raise Exception(f'Unexpected x entry')
        B = []
        for iy in range(len(a)-1):
            l = a[1+iy]
            ww = l.split()
            if len(ww) != 2:
                raise Exception(f'Unexpected data line')
            y = float(ww[0])
            if ix==0 and iy==0:
                pass
            elif ix==0 and iy==1:
                Nb = round(1/y)
            else:
                if round(Nb*y) != iy:
                    raise Exception(f'Unexpected y entry')
            B.append(int(ww[1]))
        D2.append(B)
    return Nb, D2

def covered_squares(D2, P, ix, iy):
    """
    Returns list of coordinates of square tiles that may be covered by Taylor expansion around b-lattice point ix,iy.
    """
    d2 = D2[ix][iy]
    if d2==0:
        raise Exception('Cannot add expansion as d2=0')
    pat = P[ix%2][iy%2][d2]

    ret = []
    mx = ix//2
    my = iy//2
    lx = len(pat)
    for nx in range(lx):
        jx = mx + nx - lx//2
        ly = pat[nx]
        for ny in range(ly):
            jy = ny + my - ly//2
            if jx>=0 and jy>=0:
                ret.append((jx,jy))
    return ret

def add_expansion(F, C, Q, ix, iy):
    n = len(C)
    C.append( (ix, iy) )

    qs = Q[ix][iy]

    count = 0
    for jx, jy in qs:
        if F[jx][jy] == -2:
            F[jx][jy] = n
            count += 1
        # print(f'cover {jx},{jy} by {n}')

    return count

def n_naked(qs, F):
    count = 0
    for jx, jy in qs:
        if F[jx][jy] == -2:
            count += 1
    return count

if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise Exception(f'Usage: {sys.argv[0]} <N> <file with x blocks with y tau d2 lines>')

    Ntay = int(sys.argv[1])

    # Load d2(x,y) from file.
    fname = sys.argv[2]
    Nb, D2 = read_d2_file(fname)
    d2max = max([max(line) for line in D2])
    tmax = int(sqrt(d2max)) + 1

    Rtot = 7
    Ndiv = Nb//2 # base square has edge 1/Ndiv
    Nrge = 8 # 0 <= x,y < Nrge
    Nax  = Nrge * Ndiv
    C = [] # expansion centers

    S = [[sorted_diameters(d2max, sx, sy) for sy in [0, 1]] for sx in [0, 1]]
    P = [[{d2:polyomino_pattern(d2, sx, sy) for d2 in S[sx][sy]} for sy in [0, 1]] for sx in [0, 1]]
    Q = [[covered_squares(D2, P, ix, iy) for iy in range(len(D2[ix]))] for ix in range(len(D2))]
    # print(f"42,9:\nS={S[0][1]}\nd2={D2[42][9]}\nP={P[0][1][D2[42][9]]}\nQ={sorted(Q[42][9])}\n")
    # print(f"42,8:\nS={S[0][0]}\nd2={D2[42][8]}\nP={P[0][0][D2[42][8]]}\nQ={sorted(Q[42][8])}\n")
    # print(f"41,9:\nS={S[1][1]}\nd2={D2[41][9]}\nP={P[1][1][D2[41][9]]}\nQ={sorted(Q[41][9])}\n")
    # sys.exit(1)

    # Set up field of base squares.
    # Value -1: outside domain.
    # Value -2: not yet covered by polyomino.
    F = [[-1 for jy in range(Nax)] for jx in range(Nax)]
    nF = 0
    for jx in range(Nax):
        for jy in range(Nax):
            if (jx)**2+(jy)**2 < (Rtot*Ndiv)**2:
                F[jx][jy] = -2
                nF += 1
    F[0][0] = -1 # Maclaurin expansion around z=0 is implemented separately in w_of_z.c

    # Expand around points on y axis.
    while True:
        # Search first square not yet covered
        for jy0 in range(1, Nax):
            if F[0][jy0] < 0:
                break
        if F[0][jy0] == -1: # not in domain
            break # all points inside domain are covered
        best_n = 0
        iy = None
        # print("DEBUG", 2*jy0+1, 2*Nax, 2*jy0+2*tmax+1)
        for iy1 in range(2*jy0+1, min(2*Nax, 2*jy0+2*tmax+1)):
            if F[0][iy1//2] == -1 or iy1 >= len(D2[0]):
                break
            qs = Q[0][iy1]
            if (0, jy0) in qs and len(qs) > best_n:
                best_n = n_naked(qs, F)
                iy = iy1
        if iy is None:
            break
        nF -= add_expansion(F, C, Q, 0, iy)

    # Expand around points on x axis.
    while True:
        # Search first square not yet covered
        for jx0 in range(1, Nax):
            if F[jx0][0] < 0:
                break
        if F[jx0][0] == -1: # not in domain
            break # all points inside domain are covered
        best_n = 0
        ix = None
        # print("DEBUG", 2*jx0+1, 2*Nax, 2*jx0+2*tmax+1)
        for ix1 in range(2*jx0+1, min(2*Nax, 2*jx0+2*tmax+1)):
            if F[ix1//2][0] == -1 or ix1 >= len(D2):
                break
            qs = Q[ix1][0]
            if (jx0, 0) in qs and len(qs) > best_n:
                best_n = n_naked(qs, F)
                ix = ix1
        if ix is None:
            break
        nF -= add_expansion(F, C, Q, ix, 0)

    # Sort remaining squares by distance from (-1,0).
    JJ = []
    for jx in range(1, Nax):
        for jy in range(1, Nax):
            if F[jx][jy] == -2:
                JJ.append((jx, jy))
    JJ = sorted(JJ, key=lambda jj: hypot(jj[0]+20*Ndiv, jj[1]))

    # Create expansion points inside the quadrant.
    while nF>0:
        # Search first square not yet covered
        for jx0, jy0 in JJ:
            if F[jx0][jy0] < 0:
                break
        best_n = 0
        ix = None
        for ix1 in range(1, len(D2)):
            for iy1 in range(1, len(D2[ix1])):
                if F[ix1//2][iy1//2] == -1 or D2[ix1][iy1] == 0:
                    continue
                qs = Q[ix1][iy1]
                if (jx0, jy0) in qs and len(qs) > best_n:
                    nn = n_naked(qs, F)
                    if nn > best_n:
                        ix, iy = ix1, iy1
                        best_n = nn
        if ix is None:
            break
        nF -= add_expansion(F, C, Q, ix, iy)

    print(f"needed {len(C)} subdomains")

    fname1 = "/tmp/w_taylor_centers.c"
    with open(fname1, "w") as f:
        print("// Created by %s on %s" % (" ".join(sys.argv), datetime.datetime.now().time()),
              file=f)
        print("static const int Centers[2*%i] = {" % len(C), file=f)
        for n in range(len(C)):
            ix, iy = C[n]
            print("%3i,%3i," % (ix, iy), file=f)
        print("};", file=f)
    print(f"wrote {fname1}")

    fname2 = "/tmp/w_taylor_cover.c"
    with open(fname2, "w") as f:
        print("// Created by %s on %s" % (" ".join(sys.argv), datetime.datetime.now().time()),
              file=f)
        print("static const int Cover[%i] = {" % Nax**2, file=f)
        for jx in range(Nax):
            for jy in range(Nax):
                print("%2i," % F[jx][jy], end="", file=f)
            print("", file=f)
        print("};", file=f)
    print(f"wrote {fname2}")

    fname3 = "/tmp/w_taylor_coeffs.c"
    with open(fname3, "w") as f:
        print("// Created by %s on %s" % (" ".join(sys.argv), datetime.datetime.now().time()),
              file=f)
        print("static const int NTay = %i;" % Ntay, file=f)
        print("alignas(64) static const double TaylorCoeffs[2*%i*%i] = {" % (Ntay,len(C)), file=f)
        for c in C:
            x = c[0]/(2*Ndiv)
            y = c[1]/(2*Ndiv)
            z = mpc(x, y)
            W = hp.wofz_taylor(z, Ntay)
            for w in W:
                print(" %s, %s," % (fut.double2hexstring(w.real), fut.double2hexstring(w.imag)),
                      end="", file=f)
            print(" // x=%8g y=%8g" % (x, y), file=f)
        print("};", file=f)
    print(f"wrote {fname3}")
