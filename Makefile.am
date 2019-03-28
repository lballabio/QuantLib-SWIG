
SUBDIRS = CSharp Java Python R Ruby Scala

EXTRA_DIST = \
	ChangeLog.txt \
	LICENSE.TXT \
	News.md \
	Readme.txt \
	autogen.sh


dist-hook:
	mkdir -p $(distdir)/SWIG
	cp -p ./SWIG/*.i $(distdir)/SWIG

