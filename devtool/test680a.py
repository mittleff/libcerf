#!/bin/env python

"""
Investigate algorithm 680.
"""

from mpmath import *
import sys
import hp_funcs as hp
import runtool as rt

mp.dps = 48
mp.pretty = True

def print_deviation(txt, wc, wr):
    print("%s: im %8g re %8g abs %8g" %
          (txt, abs((wc-wr).real/wr.real), abs((wc-wr).imag/wr.imag), abs((wc-wr)/wr)))
    print("new: %22.16e %22.16e" % (wc.real, wc.imag))
    print("ref: %22.16e %22.16e" % (wr.real, wr.imag))

def rho0(M, z):
    rho = mpc(0,0)
    for n in reversed(range(M)):
        rho = (n+1)/2/(z-rho)
    return rho

z = mpc(3,3)
h = mpc(0,0)
wr = hp.wofz(z)
print(wr)
M = 20

r00 = rho0(100, z)
for M in reversed(range(1,30)):
    r0 = rho0(M, z)
    w = mpc(0,1)/sqrt(pi)/(z-r0)
    print("%2i  %22.16e %22.16e  %22.16e %22.16e" %
          (M, w.real, w.imag, (w-wr).real/wr.real, (w-wr).imag/wr.imag))
