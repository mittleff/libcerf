#!/bin/env python

"""
Determine minimum of |w(z)| for z on circle with radius t around center z0.
"""

from mpmath import *
import sys
sys.path.insert(0, '../shared')
import hp_funcs as hp

mp.dps = 48
mp.pretty = True

def dw2dp(z0, t, phi):
    """
    Returns derivative d |w^2(z)| / d phi for z = z0 + t * exp(i phi).
    """
    dz = t * exp(mpc(0, phi))
    z = z0 + dz
    return 2 * re(hp.wofz(-conj(z)) * (-2*z*hp.wofz(z) + 2*mpc(0,1)/sqrt(pi)) * mpc(0,1) * dz)


def minimize_wabs(z0, t):
    """
    Determines minimum of |w(z)| for z on circle with radius t around center z0.
    """
    phi = arg(z0)
    phi = findroot(lambda phi: dw2dp(z0, t, phi), phi, solver='newton')

    dz = t * exp(mpc(0, phi))
    z = z0 + dz

    return z, fabs(hp.wofz(z)**2)

if __name__ == '__main__':
    if len(sys.argv)!=4:
        print("Usage: %s x y radius" % sys.argv[0])
        sys.exit(0)

    center = mpc(sys.argv[1], sys.argv[2])
    radius = mpf(sys.argv[3])

    zm, vm = minimize_wabs(center, radius)
    print("min at %13.4e %13.4e val %13.4e" % (zm.real, zm.imag, vm))
