#!/bin/env python

"""
Check correctness of certain functions, especially at argument values where the algorithm changes.
"""

from mpmath import *
import subprocess, sys
import hp_funcs as hp

mp.dps = 48
mp.pretty = True

### C function to be tested.

def compute_at(r):
    global run_fct_name
    mp.bps = 53
    x = mpf(r)
    fname = run_fct_name
    xs = "%22.16e" % x
    a1 = subprocess.run(['run/run_' + fname, xs], stdout=subprocess.PIPE)
    try:
        a2 = a1.stdout.decode('utf-8').split()
        a3 = [mpf(a2[0]), mpf(a2[1]), int(a2[2]), int(a2[3])]
    except:
        print("x:", x)
        print("a:", a1)
        print("Could not read back from C call")
        sys.exit(1)
    if a3[0] != mpf(xs):
        raise Exception(f"failed double-string cycle {r} -> {x} -> {xs} -> {a3[0]} ({(r-a3[0])/r})")
    mp.dps = 48
    return a3

### Bisection.

def check_at(locus, r):
    global hp_f, output_mode, worst_x, worst_relerr
    rr, f, a , n = compute_at(r)
    f2 = hp_f(rr)
    F = '%2i %3i %3i  %21.16e %21.16e  %8e %8e'
    relerr = abs(f-f2)/f2
    if relerr > worst_relerr:
        worst_x = rr
        worst_relerr = relerr
    if 't' in output_mode:
        print(F % (locus, a, n, rr, f, (f-f2)/f2, relerr))

def bisect(range_mode, r0, a0, n0, r2, a2, n2):
    global x_range
    assert(r0 < r2)
    if r2-r0 < 2e-15*(abs(r0)+abs(r2)):
        check_at(-1, r0)
        check_at(1, r2)
        return
    if range_mode == 'l':
        if r2-r0 < 1e-20 * x_range:
            check_at(-1, r0)
            check_at(1, r2)
            return
        r1 = (r0 + r2) / 2
    elif range_mode == 'p':
        r1 = sqrt(r0*r2)
    elif range_mode == 'n':
        r1 = -sqrt(r0*r2)
    wr, wi, a1, n1 = compute_at(r1)
    if (a0 != a1 or n0 != n1):
        bisect(range_mode, r0, a0, n0, r1, a1, n1)
    if (a1 != a2 or n1 != n2):
        bisect(range_mode, r1, a1, n1, r2, a2, n2)

if __name__ == '__main__':
    global output_mode, hp_f, run_fct_name, x_range, worst_x, worst_relerr
    worst_relerr = 0
    worst_x = 0

    if len(sys.argv)<7:
        print(f"Usage: {sys.argv[0]} <fct> <range_mode> <N> <from> <to> <output_mode>")
        print(f"   where <fct>         is any of i (im_wofx) x (erfcx) xm(erfcx(-x))")
        print(f"   where <range_mode>  is any of l (lin) p(pos log) n(neg log)")
        print(f"   where <output_mode> is any of t (tabulate) w(print_worst)")
        sys.exit(-1)
    if sys.argv[1] == 'i':
        run_fct_name = "imwx"
        hp_f = hp.imwx
    elif sys.argv[1] == 'x':
        run_fct_name = "erfcx"
        hp_f = hp.erfcx
    else:
        raise Exception("Invalid fct")
    range_mode = sys.argv[2]
    Ni = int(sys.argv[3])
    x_fr = float(sys.argv[4])
    x_to = float(sys.argv[5])
    output_mode = sys.argv[6]

    assert(x_fr < x_to)

    if range_mode == 'l':
        x_range = x_to-x_fr
        step = x_range/(Ni-1)
        X = [x_fr + i*step for i in range(Ni)]
    elif range_mode == 'p':
        assert(0 < x_fr)
        step = log10(x_to/x_fr)/(Ni-1)
        X = [x_fr * 10**(i*step) for i in range(Ni)]
    elif range_mode == 'n':
        assert(x_to < 0)
        step = log10(x_to/x_fr)/(Ni-1)
        X = [x_fr * 10**(i*step) for i in range(Ni)]

    r0 = None
    a0 = None
    n0 = None

    for i in range(Ni):
        r2 = X[i]
        wr, wi, a2, n2 = compute_at(r2)
        check_at(0, r2)
        if i > 0 and (a2 != a0 or n2 != n0):
            bisect(range_mode, r0, a0, n0, r2, a2, n2)
        r0 = r2
        a0 = a2
        n0 = n2
    if 'w' in output_mode:
        print("worst: at x=%22.16e relerr=%8e" % (worst_x, worst_relerr))
        rr, f, a , n = compute_at(worst_x)
        print("   run_imwx:", f)
        f2 = hp_f(worst_x, True)
        print("   highprec:", f2)
        print("   dx_rel:  %g" % ((rr-worst_x)/worst_x))
        print("   dy_rel:  %g" % ((f-f2)/f2))
        print("   algo:    ", a, n)
