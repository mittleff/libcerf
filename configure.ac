# Library libcerf:
#   Compute complex error functions, based on a new implementation of
#   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
# 
# File configure.ac:
#   Input for the build process using GNU autotools.
# 
# Copyright:
#   (C) 2012 Massachusetts Institute of Technology
#   (C) 2013 Forschungszentrum Jülich GmbH
# 
# Licence:
#   MIT Licence. See ../COPYING.
# 
# Authors:
#   Steven G. Johnson, Massachusetts Institute of Technology, 2012, core author
#   Joachim Wuttke, Forschungszentrum Jülich, 2013, package maintainer
# 
# Website:
#   http://apps.jcns.fz-juelich.de/libcerf

AC_INIT([cerf], [0:0:0], [http://apps.jcns.fz-juelich.de/libcerf])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_SRCDIR([lib/w_of_z.c]) # assure we have the right src dir
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])
# AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_SUBST(AM_CFLAGS,"-pedantic -Wall")

LT_INIT
AC_PROG_LIBTOOL
AC_PROG_CC
AC_CHECK_LIB([m], [sin], , [AC_MSG_ERROR(libm not found)])
AC_CHECK_LIB([m], [erfc], , [AC_MSG_ERROR(libm outdated: missing erfc)])
AC_HEADER_STDC
AC_CONFIG_FILES([Makefile lib/Makefile man/Makefile test/Makefile])
AC_OUTPUT
