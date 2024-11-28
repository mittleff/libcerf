#!/bin/env python

# File pro_erfcx_chebcoeffs.py:
#   Compute Chebyshev coefficients for erfcx, and write tables for use in erfcx.c.
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
mp.pretty = True

final = False # Extra checks, to be turned on in final production run

if __name__ == '__main__':

    Nout = 9

    R = fut.octavicRanges(.125, 16, 6)
    C = fut.chebcoeffs(R, Nout, hp.erfcx, final)
    fut.print_powerseries_code(R, C, Nout)
