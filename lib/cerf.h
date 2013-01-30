/* Library libcerf:
 *   compute complex error functions,
 *   along with Dawson, Faddeeva and Voigt functions
 *
 * File cerf.h:
 *   Declare exported functions.
 * 
 * Copyright:
 *   (C) 2012 Massachusetts Institute of Technology
 *   (C) 2013 Forschungszentrum Jülich GmbH
 * 
 * Licence:
 *   Permission is hereby granted, free of charge, to any person obtaining
 *   a copy of this software and associated documentation files (the
 *   "Software"), to deal in the Software without restriction, including
 *   without limitation the rights to use, copy, modify, merge, publish,
 *   distribute, sublicense, and/or sell copies of the Software, and to
 *   permit persons to whom the Software is furnished to do so, subject to
 *   the following conditions:
 * 
 *   The above copyright notice and this permission notice shall be
 *   included in all copies or substantial portions of the Software.
 * 
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *   LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *   OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 *
 * Authors:
 *   Steven G. Johnson, Massachusetts Institute of Technology, 2012, core author
 *   Joachim Wuttke, Forschungszentrum Jülich, 2013, package maintainer
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/libcerf
 *
 * Revision history:
 *   ../CHANGELOG
 *
 * More information:
 *   man 3 faddeeva
 */

#ifndef __FADDEEVA_H
#  define __FADDEEVA_H
#  undef __BEGIN_DECLS
#  undef __END_DECLS
#  ifdef __cplusplus
#    define __BEGIN_DECLS extern "C" {
#    define __END_DECLS }
#  else
#    define __BEGIN_DECLS
#    define __END_DECLS
#  endif
__BEGIN_DECLS

#include <complex.h> // C99 complex-number support

// compute w(z) = exp(-z^2) erfc(-iz) [ Faddeeva / scaled complex error func ]
extern double complex w_of_z   (double complex z,double relerr);
extern double         im_w_of_x(double x); // special case Im[w(x)] of real x

// Various functions that we can compute with the help of w(z)

// compute erfcx(z) = exp(z^2) erfc(z)
extern double complex cerfcx   (double complex z, double relerr);
extern double         faddeeva_erfcx_re(double x); // special case for real x

// compute erf(z), the error function of complex arguments
extern double complex cerf   (double complex z, double relerr);

// compute erfi(z) = -i erf(iz), the imaginary error function
extern double complex cerfi   (double complex z, double relerr);
extern double         faddeeva_erfi_re(double x); // special case for real x

// compute erfc(z) = 1 - erf(z), the complementary error function
extern double complex cerfc   (double complex z, double relerr);

// compute Dawson(z) = sqrt(pi)/2  *  exp(-z^2) * erfi(z)
extern double complex cdawson   (double complex z, double relerr);
extern double         dawson_re(double x); // special case for real x

// compute Voigt function
extern double voigt( double x );

__END_DECLS
#endif /* __FADDEEVA_H__ */
