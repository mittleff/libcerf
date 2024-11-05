#!/bin/env python

"""
Compute Taylor coefficients of w(z).
"""

from mpmath import *
import print_taylor_radius as prt

mp.dps = 48
mp.pretty = True

def per_d2(d2):
    q = lambda z: (2*prt.z2radius(z)[0])**2
    f = lambda a: mpc(0,a)
    y = findroot(lambda y: q(f(y)) - d2, [.01, 6.99], solver='bisect')
    z = f(y)
    phi0 = 0

    print(d2)
    while z.imag>=0:
        print("%8g %8g" %(z.real, z.imag))
        step = 0.1

        f = lambda a: z + step*exp(mpc(0,a))
        phi = findroot(lambda phi: q(f(phi)) - d2,
                       [phi0-pi/8, phi0+pi/8], solver='bisect')
        z = f(phi)
        phi0 = phi
    f = lambda a: mpc(a,0)
    x0 = z.real
    x = findroot(lambda x: q(f(x)) - d2, [x0-step, x0+step], solver='bisect')
    z = f(x)
    print("%8g %8g" %(z.real, z.imag))

    print()

if __name__ == '__main__':
    for d2 in [2,6,20,60]:
        per_d2(d2)
