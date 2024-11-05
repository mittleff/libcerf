#!/bin/env python

"""
Compute Taylor coefficients of w(z).
"""

from mpmath import *
import print_taylor_radius as prt

mp.dps = 48
mp.pretty = True

def per_d2(d2):
    print(d2)

    q = lambda z: (2*prt.z2radius(z)[0])**2
    f = lambda a: mpc(0,a)
    y = findroot(lambda y: q(f(y)) - d2, [.01, 2*d2], solver='bisect')
    z = f(y)
    print("%8g %8g" %(z.real, z.imag))

    z0, phi0 = z, 0
    step = min(d2/10, .05)
    while z.imag>=step:
        f = lambda a: z0 + step*exp(mpc(0,a))
        phi = findroot(lambda phi: q(f(phi)) - d2, [phi0-pi/4, phi0+pi/4], solver='bisect')
        z = f(phi)
        if abs(z) > 7:
            f = lambda a: z0 + a*exp(mpc(0,phi))
            a = findroot(lambda a: q(f(a)) - d2, [0, step], solver='bisect')
            z = f(a)
            print("%8g %8g" %(z.real, z.imag))
            print()
            return
        z = f(phi)
        print("%8g %8g" %(z.real, z.imag))
        z0, phi0 = z, phi

    try:
        f = lambda a: mpc(a,0)
        x0 = z.real
        x = findroot(lambda x: q(f(x)) - d2, [x0-5*step, x0+5*step], solver='bisect')
        z = f(x)
        print("%8g %8g" %(z.real, z.imag))
    except:
        pass
    print()

if __name__ == '__main__':
    for d2 in [.3,.5,.8,1.4,2,2.7,3.6,4.6,5.8,7.2]:
        per_d2(d2)
