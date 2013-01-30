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
