#!/bin/env python

# File pro_imwofx_chebcoeffs.py:
#   Compute Chebyshev coefficients for Im w(x), and write tables for use in im_w_of_x.c.
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
# Website:
#   http://apps.jcns.fz-juelich.de/libcerf
#
# Revision history:
#   September 2024, initial version.

from mpmath import *
import functool as fut
import random, sys

mp.dps = 48
mp.pretty = True

final = False # Extra checks, to be turned on in final production run

def hp_imwx(x, doublecheck=False):
    fz = exp(-x**2)*erfc(mpc(0, -x))
    result = fz.imag
    if doublecheck:
        # Check mpmath-computed reference value against mpmath-based brute-force integration
        r2 = 2/sqrt(pi) * exp(-x**2) * quad(lambda t : exp(t**2), [0, x])
        if abs(result-r2)/result > 1e-17:
            raise Exception(f"mpmath inaccurate")
    return fz.imag

if __name__ == '__main__':

    Nout = 8

    R = fut.octavicRanges(.5, 12., 6)
    C = fut.chebcoef(R, Nout, hp_imwx, final)

    # fut.print_cheby_coeffs(C)
    # fut.print_clenshaw_code(R, C, Nout)
    fut.print_powerseries_code(R, C, Nout)
