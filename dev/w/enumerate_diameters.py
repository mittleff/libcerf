#!/bin/env python

"""
Enumerate maximal polyominoes in contact with circumscribing circles
of increasing diameter.
"""

from math import sqrt
import sys
sys.path.insert(0, '../shared')
import functool as fut

def sorted_diameters(N):
    D2 = set()
    for j in range(1,N):
        for i in range(1,N):
            d2 = j**2 + i**2
            if d2 > N:
                continue
            D2.add(d2)
    return sorted(D2)

if __name__ == '__main__':
    fut.print_provenience()
    D2 = sorted_diameters(1450)
    for d2 in D2:
        print("%i," % d2, end="")
    print()
