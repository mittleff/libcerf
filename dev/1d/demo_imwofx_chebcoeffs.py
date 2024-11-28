#!/bin/env python

# File demo_imwofx_1cheb.py:
#   Prints Chebyshev coefficients for Im w(x) for different ranges.
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
import hp_funcs as hp

mp.dps = 48

if __name__ == '__main__':

    R = [(.5, 12, 0, 0), (.5, 1., 0, 0), (.5, .55, 0, 0), (.5, .505, 0, 0)]
    C = fut.chebcoeffs(R, 65, hp.imwx, False, 2**-51)
    fut.print_cheby_coeffs(R, C)
