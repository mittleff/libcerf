#!/bin/env python

"""
Check correctness of certain functions, especially at argument values where the algorithm changes.
"""

from mpmath import *
import sys
import hp_funcs as hp
import runtool as rt

mp.dps = 48
mp.pretty = True

if __name__ == '__main__':
    rt.worst_relerr = 0
    rt.worst_x = 0

    if len(sys.argv)<7:
        print(f"Usage: {sys.argv[0]} <fct> <range_mode> <N> <from> <to> <output_mode>")
        print(f"   where <fct>         is any of i (im_wofx) x (erfcx) xm(erfcx(-x))")
        print(f"   where <range_mode>  is any of l (lin) p(pos log) n(neg log)")
        print(f"   where <output_mode> is any of t (tabulate) w(print_worst)")
        sys.exit(-1)
    if sys.argv[1] == 'i':
        rt.run_fct_name = "imwx"
        rt.hp_f = hp.imwx
    elif sys.argv[1] == 'x':
        rt.run_fct_name = "erfcx"
        rt.hp_f = hp.erfcx
    else:
        raise Exception("Invalid fct")
    range_mode = sys.argv[2]
    Ni = int(sys.argv[3])
    x_fr = float(sys.argv[4])
    x_to = float(sys.argv[5])
    rt.output_mode = sys.argv[6]

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
        wr, wi, a2, n2 = rt.compute_at(r2)
        rt.check_at(0, r2)
        if i > 0 and (a2 != a0 or n2 != n0):
            rt.bisect(range_mode, r0, a0, n0, r2, a2, n2)
        r0 = r2
        a0 = a2
        n0 = n2
    if 'w' in rt.output_mode:
        print("rt.worst: at x=%22.16e relerr=%8e" % (rt.worst_x, rt.worst_relerr))
        rr, f, a , n = rt.compute_at(rt.worst_x)
        print("   f(x):", f)
        f2 = rt.hp_f(rt.worst_x, True)
        print("   highprec:", f2)
        print("   dx_rel:  %g" % ((rr-rt.worst_x)/rt.worst_x))
        print("   dy_rel:  %g" % ((f-f2)/f2))
        print("   algo:    ", a, n)
