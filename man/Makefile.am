# Library libcerf:
#   Compute complex error functions, based on a new implementation of
#   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
# 
# File man/Makefile.am:
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

man_MANS = faddeeva.3 voigt.3 dawson.3
auxdir = $(datadir)/man/html
aux_DATA = faddeeva.html voigt.html dawson.html

faddeeva.3 : faddeeva.pod
	pod2man -s 3 -c "libfaddeeva manual" $< > $@

voigt.3 : voigt.pod
	pod2man -s 3 -c "libfaddeeva manual" $< > $@

dawson.3 : dawson.pod
	pod2man -s 3 -c "libfaddeeva manual" $< > $@

faddeeva.html : faddeeva.pod
	pod2html --title="libfaddeeva - Faddeeva, Voigt, and Dawson functions" --noindex $< > $@

voigt.html : voigt.pod
	pod2html --title="libfaddeeva - Faddeeva, Voigt, and Dawson functions" --noindex $< > $@

dawson.html : dawson.pod
	pod2html --title="libfaddeeva - Faddeeva, Voigt, and Dawson functions" --noindex $< > $@
