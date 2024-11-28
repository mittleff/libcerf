#!/bin/env python

"""
Print coefficients exp(-a^2*n^2) for algorithm ...
"""

from mpmath import *

dps = 36
pretty = True


if __name__ == '__main__':

    a = mpf(0.5)
    a2 = a**2

    for n in range(1,100):
        y = exp(-a2*n**2)
        if y < 1e-306:
            break
        print("    %s," % (nstr(y,18,min_fixed=0)))
