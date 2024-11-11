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
    y = findroot(lambda y: f(mpc(0,y))-v, [0, 7], solver='bisect')
    z = mpc(0, y)

    ret.append((z.real, z.imag))

    z0, phi0 = z, 0
    step = .02
    while z.imag>=0:
        a = lambda phi: z0 + step*exp(mpc(0,phi))
        phi = findroot(lambda phi: f(a(phi)) - v, phi0, solver='newton')
        z = a(phi)
        if z.imag<0:
            x = findroot(lambda x: f(mpc(x,0)) - v, [0, 2*z.real], solver='bisect')
            z = mpc(x,0)
            ret.append((z.real, z.imag))
            break
        ret.append((z.real, z.imag))
        z0, phi0 = z, phi

    return ret

if __name__ == '__main__':
    for v in [.14,.3,.5,.6,.7,.75,.77,.79,.81]:
        # [.14,.165,.2,.25,.32,.4,.55,.75,1.,1.4,2.,3.,6.,15.,40.,100,250]:
        raw = per_v(v)
        out = rdp.rdp(raw, .001)

        print(v)
        for o in out:
            print("%10.5f %10.5f" % (o[0], o[1]))
        print()
