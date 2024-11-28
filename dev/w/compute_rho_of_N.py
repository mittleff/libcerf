#!/bin/env python

"""
Resolve rho**N/(1-2\rho)=epsilon for rho
"""

from mpmath import *

mp.dps = 48
mp.pretty = True

q = sqrt(5)

for N in range(8,40,4):
    f = lambda r: 2**-53*(1-2*r) - r**N
    rho = findroot(f, [.01, .49], solver='bisect')
    E = (1+(q-2)*rho-q*rho**2) / (1-(2+q)*rho+2*q*rho**2)
    print("%2i & %6.4f & %7.4g \\\\" % (N, rho, E))
