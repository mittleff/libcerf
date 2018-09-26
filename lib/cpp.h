/* Library libcerf:
 *   compute complex error functions,
 *   along with Dawson, Faddeeva and Voigt functions
 *
 * File defs.h:
 *   Define cmplx, CMPLX, NaN, for internal use, for when sources are compiled as C++ code
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

#include <cfloat>
#include <cmath>
#include <limits>
using namespace std;

// use std::numeric_limits, since 1./0. and 0./0. fail with some compilers (MS)
#define Inf numeric_limits<double>::infinity()
#define NaN numeric_limits<double>::quiet_NaN()

typedef complex<double> cmplx;

// Use C-like complex syntax, since the C syntax is more restrictive
#define cexp(z) exp(z)
#define creal(z) real(z)
#define cimag(z) imag(z)
#define cpolar(r,t) polar(r,t)

#define C(a,b) cmplx(a,b)

#define FADDEEVA(name) Faddeeva::name
#define FADDEEVA_RE(name) Faddeeva::name

// isnan/isinf were introduced in C++11
#if (__cplusplus < 201103L) && (!defined(HAVE_ISNAN) || !defined(HAVE_ISINF))
static inline bool my_isnan(double x) { return x != x; }
#  define isnan my_isnan
static inline bool my_isinf(double x) { return 1/x == 0.; }
#  define isinf my_isinf
#elif (__cplusplus >= 201103L)
// g++ gets confused between the C and C++ isnan/isinf functions
#  define isnan std::isnan
#  define isinf std::isinf
#endif

// copysign was introduced in C++11 (and is also in POSIX and C99)
#if defined(_WIN32) || defined(__WIN32__)
#  define copysign _copysign // of course MS had to be different
#elif defined(GNULIB_NAMESPACE) // we are using using gnulib <cmath>
#  define copysign GNULIB_NAMESPACE::copysign
#elif (__cplusplus < 201103L) && !defined(HAVE_COPYSIGN) && !defined(__linux__) && !(defined(__APPLE__) && defined(__MACH__)) && !defined(_AIX)
static inline double my_copysign(double x, double y) { return x<0 != y<0 ? -x : x; }
#  define copysign my_copysign
#endif

// If we are using the gnulib <cmath> (e.g. in the GNU Octave sources),
// gnulib generates a link warning if we use ::floor instead of gnulib::floor.
// This warning is completely innocuous because the only difference between
// gnulib::floor and the system ::floor (and only on ancient OSF systems)
// has to do with floor(-0), which doesn't occur in the usage below, but
// the Octave developers prefer that we silence the warning.
#ifdef GNULIB_NAMESPACE
#  define floor GNULIB_NAMESPACE::floor
#endif
