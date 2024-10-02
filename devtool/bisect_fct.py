#!/bin/env python

"""
Computes accuracy of external implementation of selected function
by comparison with an internal high-precision function.

The comparison is done for x from a linear or logarithmic grid,
and for additional points obtained by bisection around points
where the employed algorithm or number of terms has changed.

Output is either a table with one line for each x,
or a report on the worst case (the x with largest inaccuracy of f(x)).
"""

from mpmath import *
import sys
import hp_funcs as hp
import runtool as rt

mp.dps = 48
mp.pretty = True

rt.f_ext = rt.external_function1d

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
        rt.external_program = "run/run_imwx"
        rt.f_hp = hp.imwx
    elif sys.argv[1] == 'x':
        rt.external_program = "run/run_erfcx"
        rt.f_hp = hp.erfcx
    else:
        raise Exception("Invalid fct")
    rt.range_mode = sys.argv[2]
    Ni = int(sys.argv[3])
    x_fr = float(sys.argv[4])
    x_to = float(sys.argv[5])
    rt.output_mode = sys.argv[6]

    assert(x_fr < x_to)

    if rt.range_mode == 'l':
        x_range = x_to-x_fr
        step = x_range/(Ni-1)
        X = [x_fr + i*step for i in range(Ni)]
    elif rt.range_mode == 'p':
        assert(0 < x_fr)
        step = log10(x_to/x_fr)/(Ni-1)
        X = [x_fr * 10**(i*step) for i in range(Ni)]
    elif rt.range_mode == 'n':
        assert(x_to < 0)
        step = log10(x_to/x_fr)/(Ni-1)
        X = [x_fr * 10**(i*step) for i in range(Ni)]

    rt.scan_and_bisect(X)
    rt.print_conclusion()
