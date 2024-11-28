#!/bin/env python

"""
For given N, compute smallest x for which remainder R_N
of asymptotic expansion of Im w(x) falls below eps*w(x).
"""

from mpmath import *

mp.dps = 64
mp.pretty = True

eps = 2**(-53)
I = mpc(0,1)

def abs_wofx(x):
    j = mpc('0', '1')
    return abs(exp(-x**2)*erfc(-j*x))

def run(N):
    x = (gamma(N+3./2)/sqrt(pi)/eps)**(.5/(N+1)) # first order of asymptotic expansion
    # iterate using exact w(x) instead of the above estimate
    for iter in range(100):
        x = (gamma(N+3./2)/pi/eps/abs_wofx(x))**(1./(2*N+3))
    print(N, x)

if __name__ == '__main__':
    for N in range(0,50):
        run(N)
