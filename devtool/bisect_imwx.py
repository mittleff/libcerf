#!/bin/env python

"""
Check correctness of im_w_of_x, especially at argument values where the algorithm changes.
"""

import subprocess, sys
from mpmath import *

mp.dps = 48
mp.pretty = True

def dawson_kernel(t):
    return exp(t**2)

def highprecision_imwx(x, doublecheck=False):
    fz = exp(-x**2)*erfc(mpc(0, -x))
    result = fz.imag
    if doublecheck:
        # Check mpmath-computed reference value against mpmath-based brute-force integration
        r2 = 2/sqrt(pi) * exp(-x**2) * quad(dawson_kernel, [0, x])
        if abs(result-r2)/result > 1e-17:
            raise Exception(f"mpmath inaccurate")
    return fz.imag

def cerf_imwx(x):
    result = subprocess.run(['run/run_imwx', f'{x}'], stdout=subprocess.PIPE)
    a = result.stdout.decode('utf-8').split()
    return (mpf(a[0]), mpf(a[1]), int(a[2]), int(a[3]))

def compute_at(r):
    return cerf_imwx(r)

def check_at(locus, r):
    global mode, worst_x, worst_relerr
    rr, f, a , n = compute_at(r)
    f2 = highprecision_imwx(rr)
    F = '%2i %3i %3i  %12g %12g  %8e %8e'
    relerr = abs(f2-f)/f2
    if relerr > worst_relerr:
        worst_x = rr
        worst_relerr = relerr
    if 't' in mode:
        print(F % (locus, a, n, rr, f, (f2-f)/f2, relerr))

def bisect(r0, a0, n0, r2, a2, n2):
    if abs(r2-r0)<2e-15*(abs(r0)+abs(r2)):
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
    global mode, worst_x, worst_relerr
    worst_relerr = 0
    worst_x = 0

    if len(sys.argv)<5:
        print(f"Usage: {sys.argv[0]} <mode> <N> <from> <to>")
        print(f"   where <mode> is any of t (tabulate) w(print_worst)")
        sys.exit(-1)
    mode = sys.argv[1]
    Ni = int(sys.argv[2])
    x_fr = float(sys.argv[3])
    x_to = float(sys.argv[4])
    step = log10(x_to/x_fr)/(Ni-1)

    r0 = None
    a0 = None
    n0 = None

    for i in range(Ni):
        r2 = x_fr * 10**(i*step)
        wr, wi, a2, n2 = compute_at(r2)
        check_at(0, r2)
        if i > 0 and (a2 != a0 or n2 != n0):
            bisect(r0, a0, n0, r2, a2, n2)
        r0 = r2
        a0 = a2
        n0 = n2
    if 'w' in mode:
        print("worst: at x=%21.15g relerr=%8e" % (worst_x, worst_relerr))
        rr, f, a , n = compute_at(worst_x)
        print("   run_imwx:", f)
        f2 = highprecision_imwx(worst_x, True)
        print("   highprec:", f2)
        print("   dx_rel:  ", (rr-worst_x)/worst_x)
        print("   dy_rel:  ", (f-f2)/f2)
        print("   algo:    ", a, n)
