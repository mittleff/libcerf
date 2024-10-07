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
import sys

####################################################################################################
# Print functions
####################################################################################################

def print_begin_autogenerated():
    print("//--- The following code is generated by " + sys.argv[0] + "; do not edit")
    print("// clang-format off")

def print_end_autogenerated():
    print("// clang-format on")
    print("//--- End of autogenerated code")

####################################################################################################
# Conversions related to binary representation of floating-point numbers
####################################################################################################

def double2hexstring(x, nd=53):
    """
    Returns the hexagonal representation of x, rounded to nd binary digits.
    """
    if x==0:
        return "0x0"
    if x>0:
        sign = ""
    else:
        sign = "-"
        x= -x
    m,e = frexp(x)
    sm = int(nint(m*2**(nd))) # scaled mantissa
    if sm == 2**nd:
        dm = 2**(3+4*((nd-1)//4)) # displayed mantissa
        return "%s0x0.%xp%i" % (sign, dm, e+1)
    shift = 4*((nd-1)//4+1) - nd
    dm = sm*2**shift # displayed mantissa
    return "%s0x0.%0xp%i" % (sign, dm, e)

def round2(x, nd=16):
    """
    Rounds x to nd binary digits.
    """
    m,e = frexp(x)
    im = round(m * 2**nd)
    return mpf(im * 2**(e-nd))

####################################################################################################
# Computations for Chebyshev polynomials
####################################################################################################

def octavicRanges(a, b, n2e):
    """
    Divides range (a,b) into subranges, starting from octaves with boundaries at 2**m.
    Each octave is divided into 2**n2e subranges.
    Returns a list of tuples (a_subrange, b_subrange, index_of_octave, index_in_octave).
    """
    m, ea = frexp(a)
    m, eb = frexp(b)
    R = []
    for ir in range(ea-1, eb):
        for js in range(2**n2e):
            asu = (2**n2e+js) * 2**(ir-n2e)
            bsu = asu + 2.**(ir-n2e)
            if bsu>a and asu<b:
                R.append((asu, bsu, ir-ea+1, js))
    return R

def polynomial_coeffs(C):
    """
    Converts Chebyshev to power-series coefficients for faster computation in the final C code.
    Returns list P such that approximant is just sum_i P_i x^i.
    """
    N = len(C) - 1
    # Compute coefficients T_ni such that T_n(x) = sum_i T_ni x^i.
    # We also save zeros; this could be done better.
    T = [[1], [0, 1]]
    for n in range(2, N+1):
        t = [ -T[n-2][0] ]
        for i in range(1, n+1):
            if (i-n)%2==1:
                t.append(0)
            elif n>i:
                t.append(2*T[n-1][i-1]-T[n-2][i] )
            else:
                t.append(2*T[n-1][i-1] )
        T.append(t)

    # Compute coefficients of x^i, such that P_i = sum_n C_n T_ni.
    P = []
    for i in range(N+1):
        sum = 0
        for n in range(i, N+1):
            sum += C[n] * T[n][i]
        P.append(sum)

    return P

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

def check_cheb_interpolant(asu, bsu, Cs, hp_f, NT, limit):
    """
    Checks whether the interpolant with Chebyshev coefficients Cs
    agrees with the high-precision function hp_f
    within limit, for NT points within subrange (asu, bsu).
    NT should be incommensurate with len(Cs)
    """
    N = len(Cs) - 1 # polynomial order
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
        CD = [mpf("%+21.16e" % c) for c in Cs] # rounded to double precision
        ye = cheb(t, CD, N)
        r = abs((ye-yr)/yr)
        if r > limit:
            u2 = 0
            u1 = CD[N]
            msg = "%2i %+21.17f %+21.17f %+21.17f \n" % (N, CD[N], CD[N]-Cs[N], u1)
            for n in reversed(range(1,N)):
                u = 2*t*u1 - u2 + CD[n]
                msg += "%2i %+21.17f %+21.17f %+21.17f \n" % (n, Cs[n], CD[n]-Cs[n], u)
                u2 = u1
                u1 = u
            msg += "%2i %+21.17f %+21.17f %+21.17f \n" % (0, CD[0], CD[0]-Cs[0], t*u1 - u2 + CD[0])
            msg += "%46c %+21.17f \n" % (' ', yr)
            raise Exception("test failed: i=%i t=%e x=%+21.17f yr=%f err=%+21.17f relerr=%e\n%s" % (i, t, x, yr, ye-yr, r, msg))
        mp.dps = outside_dps

def chebcoeffs(R, Nout, hp_f, doublecheck, limit=2**-52):
    """
    Computes Chebyshev coefficients that interpolate the high-precision function hp_f
    for the range (a,b).
    """
    N = Nout

    C = []
    for rge in R:
        Cs = [0 for n in range(Nout+1)]
        asu, bsu, ir, js = rge

        halfrange = (bsu - asu) / 2
        center = (bsu + asu) / 2

        xx = [cos(n*pi/N) for n in range(N+1)]
        yy = [hp_f(center+halfrange*xx[n], doublecheck) for n in range(N+1)]

        for n in range(N+1):
            sum = (yy[0] + (-1)**n * yy[N]) / 2
            for j in range(1,N):
                sum += yy[j] * chebyt(n, xx[j])
            if n==0 or n==N:
                sum *= 1./N
            else:
                sum *= 2./N

            Cs[n] = sum

        check_cheb_interpolant(asu, bsu, Cs, hp_f, 316, limit)

        C.append(Cs)

    return C

####################################################################################################
# Printout for Chebyshev polynomials
####################################################################################################

def print_cheby_coeffs(R, C):
    """
    Prints Chebyshev coefficients. One line per polynomial.
    """
    nRge = len(R)
    assert(len(C) == nRge)
    for irge in range(nRge):
        Cs = C[irge]
        for coeff in Cs:
            print('%23.16e' % coeff, end=" ")
        asu, bsu, ir, js = R[irge]
        print("# subrange (%g..%g)" % (asu, bsu))

def print_clenshaw_code(R, C, Nout):
    """
    Prints C code that initializes lookup table for Clenshaw computation of Chebyshev interpolants.
    """
    nRge = len(R)
    assert(len(C) == nRge)

    print_begin_autogenerated()
    print("static const int N1cheb = %i; // number of coefficients = polynomial order + 1" % (Nout+1))
    print("static const double ChebCoeffs[%i * %i] = {" % (nRge, Nout+1))
    for irge in range(nRge):
        Cs = C[irge]
        print("   ", end="")
        for c in Cs:
            print(" %s," % double2hexstring(c), end="")
        for i in range(len(Cs), Nout+1):
            print(" 0,", end="")
        asu, bsu, ir, js = R[irge]
        print(" // x in subrange %i:%i (%g..%g)" % (ir, js, asu, bsu))
    print("};")
    print_end_autogenerated()

def print_powerseries_code(R, C, Nout):
    """
    Prints C code that initializes lookup table for Chebyshev polynomials as power series in r.
    """
    nRge = len(R)
    assert(len(C) == nRge)
    P = [polynomial_coeffs(Cs) for Cs in C]

    print_begin_autogenerated()
    koffset = 0
    if Nout + 1 == 9:
        print("alignas(64) static const double ChebCoeffs0[%i] = {" % (nRge))
        print("   ", end="")
        for irge in range(nRge):
            asu, bsu, ir, js = R[irge]
            print("     %s, // x in subrange %i:%i (%g..%g)" %
                  (double2hexstring(P[irge][0]), ir, js, asu, bsu))
        print("};")
        koffset = 1
    elif Nout + 1 == 10:
        print("alignas(64) static const double ChebCoeffs0[%i * 2] = {" % (nRge))
        print("   ", end="")
        for irge in range(nRge):
            asu, bsu, ir, js = R[irge]
            print("     %s, %s, // x in subrange %i:%i (%g..%g)" %
                  (double2hexstring(P[irge][0]), double2hexstring(P[irge][1]), ir, js, asu, bsu))
        print("};")
        koffset = 2

    print("alignas(64) static const double ChebCoeffs1[%i * 8] = {" % (nRge))
    for irge in range(nRge):
        Ps = P[irge]
        for p in Ps[koffset:Nout+1]:
            print(" %s," % double2hexstring(p), end="")
        asu, bsu, ir, js = R[irge]
        print(" // x in subrange %i:%i (%g..%g)" % (ir, js, asu, bsu))
    print("};")
    print_end_autogenerated()
