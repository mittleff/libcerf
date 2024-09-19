# File functool.py:
#   Python module providing generic tools for high-precision function evaluation.
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

####################################################################################################
# Chebyshev polynomials
####################################################################################################

def cheb(t, C, N):
    """
    Evaluates Chebyshev series at point t (between -1 and +1).
    In contrast to our final C code, we here use the Clenshaw algorithm
    [e.g. Oliver, J Inst Maths Applics 20, 379 (1977].
    """
    u2 = 0
    u1 = C[N]
    for n in reversed(range(1,N)):
        u = 2*t*u1 - u2 + C[n]
        u2 = u1
        u1 = u
    return t*u1 - u2 + C[0]

def check_cheb_interpolant(asu, bsu, C, hp_f, NT, limit):
    """
    Checks whether the interpolant with Chebyshev coefficients C
    agrees with the high-precision function hp_f
    within limit, for NT points within subrange (asu, bsu).
    NT should be incommensurate with len(C)
    """
    N = len(C) - 1 # polynomail order
    halfrange = (bsu - asu) / 2
    center = (bsu + asu) / 2
    outside_dps = mp.dps
    assert(NT>1)
    for i in range(NT):
        t = cos(i*pi/(NT-1))
        x = center+halfrange*t
        yr = hp_f(x)
        mp.dps = 16
        t = mpf(t)
        CD = [mpf("%+21.16e" % c) for c in C]
        ye = cheb(t, CD, N)
        r = abs((ye-yr)/yr)
        if r > limit:
            u2 = 0
            u1 = CD[N]
            msg = "%2i %+21.17f %+21.17f %+21.17f \n" % (N, CD[N], CD[N]-C[N], u1)
            for n in reversed(range(1,N)):
                u = 2*t*u1 - u2 + CD[n]
                msg += "%2i %+21.17f %+21.17f %+21.17f \n" % (n, C[n], CD[n]-C[n], u)
                u2 = u1
                u1 = u
            msg += "%2i %+21.17f %+21.17f %+21.17f \n" % (0, CD[0], CD[0]-C[0], t*u1 - u2 + CD[0])
            msg += "%46c %+21.17f \n" % (' ', yr)
            raise Exception("test failed: i=%i t=%e x=%+21.17f yr=%f err=%+21.17f relerr=%e\n%s" % (i, t, x, yr, ye-yr, r, msg))
        mp.dps = outside_dps
