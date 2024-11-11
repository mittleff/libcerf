#!/bin/env python

"""
Compute Taylor coefficients of w(z).
"""

from mpmath import *
import print_taylor_radius as prt
import hp_funcs as hp

mp.dps = 48
mp.pretty = True

def per_v(v):
    print(v)

    f = lambda z: abs(hp.wofz(z))
    ymin = -5
    x = findroot(lambda x: f(mpc(x,ymin))-v, [-20, -.01], solver='bisect')
    z = mpc(x, ymin)
    print("%8g %8g" %(z.real, z.imag))

    z0, phi0 = z, pi/2
    range = [0, pi]
    step = .02
    while z.imag>=ymin:
        a = lambda phi: z0 + step*exp(mpc(0,phi))
        phi = findroot(lambda phi: f(a(phi)) - v, range, solver='bisect')
        z = a(phi)
        if z.imag<ymin:
            x = findroot(lambda x: f(mpc(x,ymin)) - v, [0, 2*z.real], solver='bisect')
            z = mpc(x,ymin)
            print("%8g %8g" %(z.real, z.imag))
            break
        print("%8g %8g" %(z.real, z.imag))
        z0, phi0 = z, phi
        range = [phi0-pi/6, phi0+pi/6]
    print()

if __name__ == '__main__':
    for v in [.15,.2,.25,.32,.4,.55,.75,1.,1.4,2.,3.,6.,15.,40.,100,250]:
        per_v(v)
