# File hp_funcs.py:
#   Python module providing high-precision functions related to w(z),
#   for use in code generation for libcerf.
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

def erfcx(x, doublecheck=False):
    result = exp(x**2)*erfc(x)
    if doublecheck:
        # Check mpmath-computed reference value against mpmath-based brute-force integration
        if x<6:
            r2 = exp(x**2) - 2/sqrt(pi) * quad(lambda t : exp(x**2-t**2), [0, x])
        else:
            r2 = 2/sqrt(pi) * quad(lambda t : exp(x**2-t**2), [x, mpf('inf')])
        if abs(result-r2)/result > 1e-17:
            raise Exception(f"mpmath inaccurate")
    return result

def imwx(x, doublecheck=False):
    fz = exp(-x**2)*erfc(mpc(0, -x))
    result = fz.imag
    if doublecheck:
        # Check mpmath-computed reference value against mpmath-based brute-force integration
        r2 = 2/sqrt(pi) * exp(-x**2) * quad(lambda t : exp(t**2), [0, x])
        if abs(result-r2)/result > 1e-17:
            raise Exception(f"mpmath inaccurate")
    return result

def wofz(x, y, doublecheck=False):
    z = mpc(x,y)
    j = mpc('0', '1')
    result = exp(-z**2)*erfc(-j*z)
    if doublecheck:
        # Check mpmath-computed reference value against mpmath-based brute-force integration
        r2 = mpc(0,1)/pi*quad(lambda t: exp(-t**2)/(z-t), [-inf, +inf])
        if abs((result-r2)/result) > 1e-17:
            raise Exception(f"mpmath inaccurate")
    return result
