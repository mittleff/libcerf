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

#ifdef CERF_NO_IEEE754
#define frexp2 frexp
#else
//! Simpler replacement for frexp from math.h, assuming that value!=0.
//! Adapted from https://github.com/dioptre/newos/blob/master/lib/libm/arch/sh4/frexp.c.
inline double frexp2(double value, int* eptr)
{
    union {
	double v;
	struct {
            unsigned u_mant2 : 32;
            unsigned u_mant1 : 20;
            unsigned   u_exp : 11;
	    unsigned  u_sign :  1;
	} s;
    } u;

    u.v = value;
    *eptr = u.s.u_exp - 1022;
    u.s.u_exp = 1022;
    return u.v ;
}
#endif
