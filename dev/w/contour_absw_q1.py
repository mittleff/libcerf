#!/bin/env python

"""
Contour of |w(z)|=const in 1st quadrant.
"""

from mpmath import *
import print_taylor_radius as prt
import hp_funcs as hp
import rdp

mp.dps = 48
mp.pretty = True

def per_v(v):
    ret = []

    f = lambda z: abs(hp.wofz(z))
    y = findroot(lambda y: f(mpc(0,y))-v, [0, 10], solver='bisect')
    z = mpc(0, y)

    ret.append((z.real, z.imag))

    z0, phi0 = z, 0
    step = .02
    while True:
        a = lambda phi: z0 + step*exp(mpc(0,phi))
        phi = findroot(lambda phi: f(a(phi)) - v, phi0, solver='newton')
        z = a(phi)
        ymin = -1
        if z.imag<ymin:
            x = findroot(lambda x: f(mpc(x,ymin)) - v, z.real, solver='newton')
            z = mpc(x,ymin)
            ret.append((z.real, z.imag))
            break
        ret.append((z.real, z.imag))
        z0, phi0 = z, phi

    return ret

if __name__ == '__main__':
    for v in [.065,.075,.088,.105,.132,.177,.247,.39,.67]:
        print(v)
        raw = per_v(v)
        out = rdp.rdp(raw, .001)
        for o in out:
            print("%10.5f %10.5f" % (o[0], o[1]))
        print()
