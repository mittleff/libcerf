#!/bin/env python

"""
Enumerate maximal polyominoes in contact with circumscribing circles
of increasing diameter.
"""

from math import sqrt

N=13
Raw = []
for j in range(1,N):
    for i in range(j,N):
        d2 = j**2 + i**2
        if d2 > N**2:
            continue
        Raw.append((i, j, d2))
Tup = sorted(Raw, key=lambda t: t[2])

Out = []
for n in range(len(Tup)):
    i,j,d2 = Tup[n]
    rows = []
    if j%2 == 0:
        for jj in range(1, int(sqrt(d2)/2+1)):
            dx = sqrt(d2/4 - jj**2)
            if (i%2 == 0):
                if dx > 0:
                    row = 2 * int(dx)
                    rows.insert(0, row)
                    rows.append(row)
            else:
                if dx >= 0.5:
                    row = 1 + 2 * int(dx - 0.5)
                    rows.insert(0, row)
                    rows.append(row)
    else:
        for jj2 in range(1, int(sqrt(d2))+1, 2):
            dx = sqrt(d2 - jj2**2)/2
            if (i%2 == 0):
                if dx >= 1:
                    row = 2 * int(dx)
                    if jj2 > 1:
                        rows.insert(0, row)
                    rows.append(row)
            else:
                if dx >= 0.5:
                    row = 1 + 2 * int(dx-0.5)
                    if jj2 > 1:
                        rows.insert(0, row)
                    rows.append(row)
    area = 0
    for r in rows:
        area += r

#    print("%2i %2i %2i %3i %3i %s" % (n, i, j, d2, area, rows))
    print("%2i [%s] P" % (d2, ' '.join([str(r) for r in rows])))
