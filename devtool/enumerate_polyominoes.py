#!/bin/env python

"""
Enumerate maximal polyominoes in contact with circumscribing circles
of increasing diameter.
"""

from math import sqrt
import functool as fut

def sorted_polyominoes(N):
    Out = []
    for s in range(1,N):
        for t in range(s,N):
            d2 = s**2 + t**2
            if d2 > N**2:
                continue

            rows = []
            nj2 = int(sqrt(d2)) + 1
            for k in range(2-s%2, nj2, 2):
                dx = sqrt(d2 - k**2)/2 - (t%2)/2
                row = 2 * int(dx) + (t%2)
                if row > 0:
                    if k > s%2:
                        rows.insert(0, row)
                    rows.append(row)
            area = 0
            for n in rows:
                area += n

            Out.append((s, t, d2, area, rows))

    return sorted(Out, key=lambda e: e[2])

if __name__ == '__main__':
    fut.print_provenience()
    P = sorted_polyominoes(10)
    for p in P:
        i, j, d2, area, rows = p
        #    print("%2i %2i %2i %3i %3i %s" % (n, i, j, d2, area, rows))
        print("%3i [%s] P" % (d2, ' '.join([str(n) for n in rows])))
    print()
