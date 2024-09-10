#!/bin/env python

"""
Check correctness of im_w_of_x, especially at argument values where the algorithm changes.
"""

import subprocess
from mpmath import *

mp.dps = 48
mp.pretty = True

def highprecision_imwx(x):
    z = mpc(x,0)
    j = mpc('0', '1')
    fz = exp(-z**2)*erfc(-j*z)
    return fz.imag

def cerf_imwx(x):
    result = subprocess.run(['run/run_imwx', f'{x}'], stdout=subprocess.PIPE)
    a = result.stdout.decode('utf-8').split()
    return (mpf(a[0]), mpf(a[1]), int(a[2]), int(a[3]))

def compute_at(r):
    return cerf_imwx(r)

def check_at(locus, r):
    rr, f, a , n = compute_at(r)
    f2 = highprecision_imwx(r)
    F = '%2i %3i %3i  %12g %12g  %8e'
    #if abs(wr2-wr) > 3e-15 * wr2:
    print(F % (locus, a, n, r, f, abs(f2-f)/f2))

def bisect(r0, a0, n0, r2, a2, n2):
    if abs(r2-r0)<2e-15*(abs(r0)+abs(r2)):
        # print(f'jump {a0}/{n0} -> {a2}/{n2} after {r0}')
        check_at(-1, r0)
        check_at(1, r2)
        return
    r1 = sqrt(r0*r2)
    wr, wi, a1, n1 = compute_at(r1)
    if (a0 != a1 or n0 != n1):
        bisect(r0, a0, n0, r1, a1, n1)
    if (a1 != a2 or n1 != n2):
        bisect(r1, a1, n1, r2, a2, n2)

if __name__ == '__main__':
    global s
    r0 = None
    a0 = None
    n0 = None
    Ni = 13241
    for i in range(Ni):
        r2 = mpf('10')**((i-Ni/2.)*2/Ni*8)
        wr, wi, a2, n2 = compute_at(r2)
        check_at(0, r2)
        if i > 0 and (a2 != a0 or n2 != n0):
            bisect(r0, a0, n0, r2, a2, n2)
        r0 = r2
        a0 = a2
        n0 = n2
