#!/bin/env python

"""
Resolve rho**N/(1-2\rho)=epsilon for rho
"""

from mpmath import *

mp.dps = 48
mp.pretty = True


for N in range(8,40,4):
    f = lambda r: 2**-53*(1-2*r) - r**N
    rho = findroot(f, [.01, .49], solver='bisect')
    E = 1/(1-rho)/(1-2*rho)
    print(N, rho, E)
