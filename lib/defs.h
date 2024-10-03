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
