# Library libcerf:
#   Compute complex error functions, based on a new implementation of
#   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
# 
# File test/Makefile.am:
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

AM_CFLAGS = -I$(srcdir)/../lib -Wall -std=c99
bin_PROGRAMS = test_libcerf_1
test_libcerf_1_SOURCES = test1.c
AM_LDFLAGS = -lm -L../lib -lcerf
