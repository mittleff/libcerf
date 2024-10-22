#!/bin/env python

"""
Compute Taylor coefficients of w(z).
"""

from mpmath import *
import sys
import hp_funcs as hp

mp.dps = 48
mp.pretty = True

def forward(z):
    W = []
    W.append(hp.wofz(z.real, z.imag))
    W.append(-2*z*W[0] + mpc(0,2)/sqrt(pi))
    for k in range(2,20):
        W.append(-2*(z*W[k-1]+(k-1)*W[k-2]))
    return W

def backward(z, M):
    w = z
    W = [0 for k in range(20)]
    for k in reversed(range(M)):
        w = z - (k+1)/2/w
        if k<20:
            W[k] = mpc(0,1)/sqrt(pi)/w
            print("bw %2i %22.15e %22.15e" % (k, W[k].real, W[k].imag))
    return W

if __name__ == '__main__':
    z = mpc(4.5, 3.5)
    Wa = forward(z)
    Wb = backward(z, 99)
    for k in range(0,20):
        print("%2i %22.15e %22.15e  %22.15e %22.15e" %
              (k, Wa[k].real, Wb[k].real, Wa[k].imag, Wb[k].imag))
