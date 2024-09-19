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
            if bsu>=a and asu<=b:
                R.append((asu, bsu, ir-ea+1, js))
    return R

def polynomial_coefs(C):
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

def chebcoef(R, Nout, hp_f, doublecheck):
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

        for m in range(N+1):
            sum = (yy[0] + (-1)**m * yy[N]) / 2
            for n in range(1,N):
                sum += yy[n] * chebyt(m, xx[n])
            if m==0 or m==N+1:
                sum *= 1./N
            else:
                sum *= 2./N

            #if ((m > 1 and not final and not tuning) or m>Nout) and abs(sum) < 5e-17*abs(C[s][0]):
            #    break
            #if m > Nout:
            #    raise Exception(f"N={Nout} exceeded")
            Cs[m] = sum

        # print('%3i %8e %8e %+8e %2i' % (s, asu, bsu, C[s][m-1]/C[s][0], m-1))
        #for mm in range(m, Nout+1):
        #    del C[s][-1]

        check_cheb_interpolant(asu, bsu, Cs, hp_f, 316, 2**-52)

        C.append(Cs)

    return C

####################################################################################################
# Printout for Chebyshev polynomials
####################################################################################################

def print_cheby_coeffs(C):
    for Cs in C:
        for coeff in Cs:
            print('%23.16e' % coeff, end=", ")
        print()
    print()

def print_clenshaw_code(R, C, Nout):
    """
    Prints C code that initializes lookup table for Clenshaw computation of Chebyshev interpolants.
    """
    nRge = len(R)
    assert(len(C) == nRge)

    print("//--- The following code is generated by " + sys.argv[0])
    print("// clang-format off")
    print("static const int N1cheb = %i; // number of coefficients = polynomial order + 1" % (Nout+1))
    print("static const double ChebCoeffs[%i * %i] = {" % (nRge, Nout+1))
    for irge in range(nRge):
        Cs = C[irge]
        print("   ", end="")
        for c in Cs:
            print(" %+22.16e," % c, end="")
        for i in range(len(Cs), Nout+1):
            print(" 0,", end="")
        asu, bsu, ir, js = R[irge]
        print(" // x in subrange %i:%i (%g..%g)" % (ir, js, asu, bsu))
    print("};")
    print("// clang-format on")
    print("//--- End of autogenerated code")

def print_powerseries_code(R, C, Nout):
    """
    Prints C code that initializes lookup table for Chebyshev polynomials as power series in r.
    """
    nRge = len(R)
    assert(len(C) == nRge)
    P = [polynomial_coefs(Cs) for Cs in C]

    print("//--- The following code is generated by " + sys.argv[0])
    print("// clang-format off")
    # print("static const int N1cheb = %i; // number of coefficients = polynomial order + 1" % (Nout+1))
    N2 = Nout + 1
    if N2 == 9:
        print("static const double ChebCoeffs8[%i] ALIGN(64) = {" % (nRge))
        print("   ", end="")
        for irge in range(nRge):
            asu, bsu, ir, js = R[irge]
            p8 = P[irge][8]
            print("     %+22.16e, // x in subrange %i:%i (%g..%g)" % (p8, ir, js, asu, bsu))
        print("};")
        N2 -= 1
    elif N2 == 10:
        print("static const double ChebCoeffs89[%i * 2] ALIGN(64) = {" % (nRge))
        print("   ", end="")
        for irge in range(nRge):
            asu, bsu, ir, js = R[irge]
            print("     %+22.16e, %+22.16e, // x in subrange %i:%i (%g..%g)" %
                  (P[irge][8], P[irge][9], ir, js, asu, bsu))
        print("};")
        N2 -= 2

    print("static const double ChebCoeffs[%i * %i] ALIGN(64) = {" % (nRge, N2))
    for irge in range(nRge):
        Cs = C[irge]
        Ps = P[irge]
        for p in Ps[0:N2]:
            print(" %+22.16e," % p, end="")
        asu, bsu, ir, js = R[irge]
        print(" // x in subrange %i:%i (%g..%g)" % (ir, js, asu, bsu))
    print("};")
    print("// clang-format on")
    print("//--- End of autogenerated code")
