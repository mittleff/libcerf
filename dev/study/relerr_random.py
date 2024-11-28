#!/bin/env python

"""
Computes relative error of an ad-hoc function for random input
"""

from mpmath import *
import sys

mp.bps = 106
mp.pretty = True

if __name__ == '__main__':
    wc_relerr = 0
    n = 0
    while True:
        n += 1

        x = 1+rand()
        y = 1+rand()
        z = mpc(x, y)
        v = 1/z

        xf = float(x)
        yf = float(y)
        zf = complex(xf, yf)
        vf = 1/zf

        re = abs((v-vf)/v)
        if re > wc_relerr:
            wc_relerr = re
            print("%8i %21.16g %21.16g %12g" % (n, x, y, re/2**-53))
