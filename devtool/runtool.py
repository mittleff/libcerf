# File runtool.py:
#   This Python module provides generic tools for
#   - executing a C program that runs a numeric computation
#   - bisection to find arguments hwere the numeric algorithm has changed
#
# Copyright:
#   (C) 2024 Forschungszentrum Jülich GmbH
#
# Licence:
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the
#   "Software"), to deal in the Software without restriction, including
#   without limitation the rights to use, copy, modify, merge, publish,
#   distribute, sublicense, and/or sell copies of the Software, and to
#   permit persons to whom the Software is furnished to do so, subject to
#   the following conditions:
#
#   The above copyright notice and this permission notice shall be
#   included in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#   LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#   OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
#   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# Author:
#   Joachim Wuttke, Forschungszentrum Jülich, 2024
#
# Revision history:
#   Initial version published in libcerf/devtool.

from mpmath import *
import subprocess, sys

this = sys.modules[__name__]

### Run C function to be tested.

def external_function1d(x):
    """
    Evaluates a one-dimensional real function f(x) by calling an external program.
    Rounds r to 53 binary digits to avoid loss of precision in external program.
    Returns tuple (x, f(x), algorithm_number, number_of_terms).
    """
    mp.bps = 53
    xhp = mpf(x)
    xs = "%22.16e" % xhp
    a1 = subprocess.run([this.external_program, xs], stdout=subprocess.PIPE)
    try:
        a2 = a1.stdout.decode('utf-8').split()
        a3 = [mpf(a2[0]), mpf(a2[1]), int(a2[2]), int(a2[3])]
    except:
        print("x:", x)
        print("a:", a1)
        print("Could not read back from C call")
        sys.exit(1)
    if a3[0] != mpf(xs):
        raise Exception(f"failed double-string cycle {x} -> {xhp} -> {xs} -> {a3[0]} ({(r-a3[0])/r})")
    mp.dps = 48
    return a3

### Bisection.

def check_at(locus, r):
    rr, f, a , n = this.f_ext(r)
    f2 = this.f_hp(rr)
    F = '%2i %3i %3i  %21.16e %21.16e  %8e %8e'
    relerr = abs(f-f2)/f2
    if relerr > this.worst_relerr:
        this.worst_x = rr
        this.worst_relerr = relerr
    if 't' in this.output_mode:
        print(F % (locus, a, n, rr, f, (f-f2)/f2, relerr))

def bisect(range_mode, r0, a0, n0, r2, a2, n2):
    assert(r0 < r2)
    if r2-r0 < 2e-15*(abs(r0)+abs(r2)):
        check_at(-1, r0)
        check_at(1, r2)
        return
    if range_mode == 'l':
        if r2-r0 < 1e-20 * this.x_range:
            check_at(-1, r0)
            check_at(1, r2)
            return
        r1 = (r0 + r2) / 2
    elif range_mode == 'p':
        r1 = sqrt(r0*r2)
    elif range_mode == 'n':
        r1 = -sqrt(r0*r2)
    wr, wi, a1, n1 = this.f_ext(r1)
    if (a0 != a1 or n0 != n1):
        bisect(range_mode, r0, a0, n0, r1, a1, n1)
    if (a1 != a2 or n1 != n2):
        bisect(range_mode, r1, a1, n1, r2, a2, n2)

def scan_and_bisect(X):
    """
    Computes relative accuracy for all x in X.
    Furthermore, does a bisection whenever algorithm or number of terms has changed.
    """
    r0 = None
    a0 = None
    n0 = None

    for i in range(len(X)):
        r2 = X[i]
        wr, wi, a2, n2 = f_ext(r2)
        check_at(0, r2)
        if i > 0 and (a2 != a0 or n2 != n0):
            bisect(range_mode, r0, a0, n0, r2, a2, n2)
        r0 = r2
        a0 = a2
        n0 = n2

# Reporting.

def print_conclusion():
    if 'w' in this.output_mode:
        print("this.worst: at x=%22.16e relerr=%8e" % (this.worst_x, this.worst_relerr))
        rr, f, a , n = this.this.f_ext(this.worst_x)
        print("   f(x):", f)
        f2 = this.f_hp(this.worst_x, True)
        print("   highprec:", f2)
        print("   dx_rel:  %g" % ((rr-this.worst_x)/this.worst_x))
        print("   dy_rel:  %g" % ((f-f2)/f2))
        print("   algo:    ", a, n)
