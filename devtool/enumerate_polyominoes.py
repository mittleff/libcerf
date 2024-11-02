#!/bin/env python

"""
Enumerate maximal polyominoes in contact with circumscribing circles
of increasing diameter.
"""

from math import sqrt
import functool as fut

def sorted_polyominoes(N):
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
            nj2 = 2*int(sqrt(d2)/2+1)
        else:
            nj2 = int(sqrt(d2))+1
        for jj2 in range(2-j%2, nj2, 2):
            dx = sqrt(d2 - jj2**2)/2 - (i%2)/2
            row = 2 * int(dx) + (i%2)
            if row > 0:
                if jj2 > j%2:
                    rows.insert(0, row)
                rows.append(row)
        area = 0
        for r in rows:
            area += r

        Out.append((n, i, j, d2, area, rows))

    return Out

if __name__ == '__main__':
    fut.print_provenience()
    P = sorted_polyominoes(10)
    for p in P:
        n, i, j, d2, area, rows = p
        #    print("%2i %2i %2i %3i %3i %s" % (n, i, j, d2, area, rows))
        print("%3i [%s] P" % (d2, ' '.join([str(r) for r in rows])))
    print()
