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
        for i in range(1,N):
            d2 = j**2 + i**2
            if d2 > N**2:
                continue
            Raw.append((i, j, d2))
    Tup = sorted(Raw, key=lambda t: t[2])
    return Tup

if __name__ == '__main__':
    fut.print_provenience()
    P = sorted_polyominoes(50)
    for p in P:
        i, j, d2 = p
        # print("%2i %2i %3i" % (i, j, d2))
        print("%i," % j, end="")
    print()
