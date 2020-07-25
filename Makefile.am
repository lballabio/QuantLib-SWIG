
SUBDIRS = CSharp Java Python R Scala

EXTRA_DIST = \
	ChangeLog.txt \
	LICENSE.TXT \
	News.md \
	README.md \
	autogen.sh


dist-hook:
	mkdir -p $(distdir)/SWIG
	cp -p ./SWIG/*.i $(distdir)/SWIG
	cp -p ./requirements.txt $(distdir)/Python/examples/
