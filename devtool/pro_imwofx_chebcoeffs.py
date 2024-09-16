#!/bin/env python

"""
Compute Chebyshev coefficients for im_w_of_x.
Used to generate code sections in im_w_of_z.c.
"""

from mpmath import *
import random, sys

mp.dps = 48
mp.pretty = True

tuning = True
final = False # Extra checks, to be turned on in final production run

def dawson_kernel(t):
    return exp(t**2)

def highprecision_imwx(x):
    fz = exp(-x**2)*erfc(mpc(0, -x))
    result = fz.imag
    if final:
        # Check mpmath-computed reference value against mpmath-based brute-force integration
        r2 = 2/sqrt(pi) * exp(-x**2) * quad(dawson_kernel, [0, x])
        if abs(result-r2)/result > 1e-17:
            raise Exception(f"mpmath inaccurate")
    return fz.imag

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

def test(asu, bsu, C):
    """
    Checks our Chebyshev interpolant against the target function,
    for lots of points within given subrange.
    """
    N = len(C) - 1
    halfrange = (bsu - asu) / 2
    center = (bsu + asu) / 2
    NT = 316 # NT+1 should be incommensurate with N+1
    for i in range(NT+1):
        t = cos(i*pi/NT)
        x = center+halfrange*t
        yr = highprecision_imwx(x)
        mp.dps = 16
        t = mpf(t)
        CD = [mpf("%+21.16e" % c) for c in C]
        ye = cheb(t, CD, N)
        r = abs((ye-yr)/yr)
        if r > 2.24e-16:
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
        mp.dps = 48

def subrange(a, b, S, s):
    asu = ((S-s)*a + s*b) / S
    bsu = ((S-s-1)*a + (s+1)*b) / S
    return (asu, bsu)

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

def chebcoef(R, Nout):
    """
    Computes Chebyshev coefficients that interpolate function highprecision_imwx in the range (a,b).
    """
    N = Nout

    C = []
    for rge in R:
        Cs = [0 for n in range(Nout+1)]
        asu, bsu, ir, js = rge

        halfrange = (bsu - asu) / 2
        center = (bsu + asu) / 2

        xx = [cos(n*pi/N) for n in range(N+1)]
        yy = [highprecision_imwx(center+halfrange*xx[n]) for n in range(N+1)]

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

        test(asu, bsu, Cs)

        C.append(Cs)

    return C

if __name__ == '__main__':

    Nout = 8

    n2e = 6
    R = []
    for ir in range(5):
        for js in range(2**n2e):
            a = (2**n2e+js) * 2**(ir-n2e-1)
            b = a + 2.**(ir-n2e-1)
            if b>=.5 and a<=12.:
                R.append((a,b, ir, js))
    C = chebcoef(R, Nout)

    # print_cheby_coeffs(C)
    # print_clenshaw_code(R, C, Nout)
    print_powerseries_code(R, C, Nout)
