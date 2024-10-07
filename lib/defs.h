/* Library libcerf:
 *   compute complex error functions,
 *   along with Dawson, Faddeeva and Voigt functions
 *
 * File defs.h:
 *   Language-dependent includes.
 *
 * Copyright:
 *   (C) 2012 Massachusetts Institute of Technology
 *   (C) 2013 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   MIT Licence.
 *   See ../COPYING
 *
 * Authors:
 *   Steven G. Johnson, Massachusetts Institute of Technology, 2012, core author
 *   Joachim Wuttke, Forschungszentrum Jülich, 2013, package maintainer
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/libcerf
 */


#ifdef __cplusplus
#include "cpp.h"
#else
#include "c.h"
#endif

#ifndef CERF_AS_CPP
#define alignas _Alignas // C23 will do this for us
#endif

#ifdef CERF_INTROSPECT
#define SET_INFO(a,n) cerf_algorithm = a; cerf_nofterms = n
#define SET_ALGO(a) cerf_algorithm = a
#define SET_NTER(n) cerf_nofterms = n
#else
#define SET_INFO(a,n)
#define SET_ALGO(a)
#define SET_NTER(n)
#endif

#ifdef CERF_NO_IEEE754 // This flag can be set via CMake option -DCERF_IEEE754=OFF
// Fall back to frexp from math.h. To be used for non-standard processor architectures
// for which our accelerated function frexp2 does not work.
#define frexp2 frexp

#else
//! Simpler replacement for frexp from math.h, assuming that 0 < value < inf.
//!
//! Adapted from https://github.com/dioptre/newos/blob/master/lib/libm/arch/sh4/frexp.c.
//!
//! This function assumes that integer and floating-point numbers have the same endianness.
//! This seems to be case "on modern standard computers" that implement IEEE 754 [1].
//! Otherwise set flag CERF_NO_IEEE754 to fall back to frexp.
//!
//! [1] https://en.wikipedia.org/wiki/Endianness#Floating_point of 07oct24
inline double frexp2(double value, int* eptr)
{
    union {
	double v;
	struct {
	    // mantissa split in two parts lest we get segfault under MSVC
            unsigned mantissa2 : 32;
            unsigned mantissa1 : 20;
            unsigned exponent : 11;
	    unsigned sign :  1;
	} s;
    } u;

    u.v = value;
    *eptr = u.s.exponent - 1022;
    u.s.exponent = 1022;
    return u.v ;
}
#endif
