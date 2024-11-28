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
mp.dps = 48
mp.pretty = True

#RadialData = (
# Created by devtool/print_taylor_radius.py on 21:01:25.585837
#(  0,   0,  17,   2,   0.2417), # [2, 2]

def neighbors(j, i, kP, Rows):
    rows = Rows[kP]
    result = []
    nj = len(rows)
    assert(j%2 == nj%2)
    assert(i%2 == rows[0]%2)
    j0 = j//2
    for jj in range(nj):
        jcorner = j0 + jj - nj//2
        if jcorner<0 or jcorner>=Nax:
            continue
        ni = rows[jj]
        i0 = i//2
        for ii in range(ni):
            icorner = i0 + ii - ni//2
            if icorner<0 or icorner>=Nax:
                continue
            result.append((jcorner, icorner))
    # print("neighbors(%3i %3i) kP=%3i nj=%3i %s -> %s (size %i)" % (j, i, kP, nj, str(rows), str(result), len(result)))
    return result

def unclaimed(F):
    result = 0
    for row in F:
        for val in row:
            if val == -2:
                result += 1
    return result

def adjust_d2(D2, ix, iy, S):
    """
    Reduces D2 entry to value with correct parities.
    """
    d2 = D2[ix][iy]
    if not d2 in S:
        for s in reversed(S):
            if s <= d2:
                D2[ix][iy] = s
                return
        D2[ix][iy] = 0

def sorted_diameters(N, sx, sy):
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
    Read file that provides d2(x,y).
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
            if len(ww) != 3:
                raise Exception(f'Unexpected data line')
            y = float(ww[0])
            if ix==0 and iy==0:
                pass
            elif ix==0 and iy==1:
                Nb = round(1/y)
            else:
                if round(Nb*y) != iy:
                    raise Exception(f'Unexpected y entry')
            B.append(int(ww[2]))
        D2.append(B)
    return Nb, D2

def covered_squares(F, D2, P, ix, iy):
    d2 = D2[ix][iy]
    if d2==0:
        raise Exception('Cannot add expansion as d2=0')
    pat = P[ix%2][iy%2][d2]

    ret = []
    mx = ix//2
    my = iy//2
    ly = len(pat)
    for ny in range(ly):
        jy = my + ny - ly//2
        lx = pat[ny]
        for nx in range(lx):
            jx = nx + mx - lx//2
            if jx>=0 and jy>=0 and F[jx][jy] == -2:
                ret.append((jx,jy))
    return ret

def add_expansion(F, C, D2, P, ix, iy):
    n = len(C)
    C.append( (ix, iy) )

    qs = covered_squares(F, D2, P, ix, iy)

    for jx, jy in qs:
        assert(F[jx][jy] == -2)
        F[jx][jy] = n
        # print(f'cover {jx},{jy} by {n}')

    return len(qs)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise Exception(f'Usage: {sys.argv[0]} <file with x blocks with y tau d2 lines>')

    # Load d2(x,y) from file.
    fname = sys.argv[1]
    Nb, D2 = read_d2_file(fname)
    d2max = max([max(line) for line in D2])
    tmax = int(sqrt(d2max)) + 1

    # Decrease d2 so that it matches the parity of the lattice point.
    S = [[sorted_diameters(d2max, sx, sy) for sy in [0, 1]] for sx in [0, 1]]
    for ix in range(len(D2)):
        B = D2[ix]
        for iy in range(len(B)):
            adjust_d2(D2, ix, iy, S[ix%2][iy%2])

    # Precompute patterns P[sx][sy][d2].
    P = [[{d2:polyomino_pattern(d2, sx, sy) for d2 in S[sx][sy]} for sy in [0, 1]] for sx in [0, 1]]

    Rtot = 7
    Ndiv = Nb//2 # base square has edge 1/Ndiv
    Nrge = 8 # 0 <= x,y < Nrge
    Nax  = Nrge * Ndiv
    C = [] # expansion centers

    # Set up field of base squares.
    # Value -1: outside domain.
    # Value -2: not yet covered by polyomino.
    F = [[-1 for jy in range(Nax)] for jx in range(Nax)]
    nF = 0
    for jx in range(Nax):
        for jy in range(Nax):
            if jx**2+jy**2 <= (Rtot*Ndiv)**2:
                F[jx][jy] = -2
                nF += 1

    fut.print_provenience()

    # Expand around (0,0).
    nF -= add_expansion(F, C, D2, P, 0, 0)

    # Expand around points on y axis.
    while True:
        # Search first point not yet covered
        for jy0 in range(1, Nax):
            if F[0][jy0] < 0:
                break
        if F[0][jy0] == -1: # not in domain
            break # all points inside domain are covered
        options = []
        best_n = 0
        iy = None
        # print("DEBUG", 2*jy0+1, 2*Nax, 2*jy0+2*tmax+1)
        for iy1 in range(2*jy0+1, min(2*Nax, 2*jy0+2*tmax+1)):
            if F[0][iy1//2] == -1 or iy1 >= len(D2[0]):
                break
            qs = covered_squares(F, D2, P, 0, iy1)
            if (0, jy0) in qs and len(qs) > best_n:
                best_n = len(qs)
                iy = iy1
        if iy is None:
            break
        nF -= add_expansion(F, C, D2, P, 0, iy)

    # Expand around points on x axis.
    while True:
        # Search first point not yet covered
        for jx0 in range(1, Nax):
            if F[jx0][0] < 0:
                break
        if F[jx0][0] == -1: # not in domain
            break # all points inside domain are covered
        options = []
        best_n = 0
        ix = None
        # print("DEBUG", 2*jx0+1, 2*Nax, 2*jx0+2*tmax+1)
        for ix1 in range(2*jx0+1, min(2*Nax, 2*jx0+2*tmax+1)):
            if F[ix1//2][0] == -1 or ix1 >= len(D2):
                break
            qs = covered_squares(F, D2, P, ix1, 0)
            if (jx0, 0) in qs and len(qs) > best_n:
                best_n = len(qs)
                ix = ix1
        if ix is None:
            break
        nF -= add_expansion(F, C, D2, P, ix, 0)

#     while unclaimed(F)>0:
#         max_improve = 0
#         i_max, j_max = -1, -1
#         for j in range(2*Nax):
#             for i in range(2*Nax):
#                 kP = RadialIndex[j][i]
#                 Neighbors = neighbors(j, i, kP, Rows)
#                 count_improve = 0
#                 for jj, ii in Neighbors:
#                     if F[jj][ii] == -2:
#                         count_improve += 1
#                 if count_improve > max_improve:
#                     max_improve = count_improve
#                     j_max = j
#                     i_max = i
#
#         kP = RadialIndex[j_max][i_max]
#         Neighbors = neighbors(j_max, i_max, kP, Rows)
#         for jj, ii in Neighbors:
#             if F[jj][ii] == -2:
#                 F[jj][ii] = len(C)
#         C.append((i_max,j_max))
#         print("# Improve by %3i Place %3i at %3i %3i, kP=%3i neighbors #=%3i" %
#               (max_improve, len(C), j_max, i_max, RadialIndex[j][i], len(Neighbors)))
#
    print("# coord -> tile")
    for jx in range(Nax):
        for jy in range(Nax):
            print("%i," % F[jx][jy], end="")
        print()
    print("# tiles")
    for c in C:
        print("%3i, %3i," % c)
#
#     print("# Taylor coefficients")
#     fut.print_begin_autogenerated()
#     print("alignas(64) static const double TaylorCoeffs[%i] = {" % (2*24*len(C)))
#     for c in C:
#         x = c[0]/(2*Ndiv)
#         y = c[1]/(2*Ndiv)
#         z = mpc(x, y)
#         W = hp.wofz_taylor(z, 24)
#         for w in W:
#             print(" %s, %s," % (fut.double2hexstring(w.real), fut.double2hexstring(w.imag)), end="")
#         print(" // x=%8g y=%8g" % (x, y))
#     print("};")
#     fut.print_end_autogenerated()
