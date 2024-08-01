
CLEANFILES = src/QuantLib/quantlib_wrap.cpp src/QuantLib/QuantLib.py .build-stamp

BUILT_SOURCES = src/QuantLib/quantlib_wrap.cpp src/QuantLib/QuantLib.py

if HAVE_PYTHON
if BUILD_PYTHON

all-local: .build-stamp

.build-stamp: src/QuantLib/quantlib_wrap.cpp src/QuantLib/QuantLib.py setup.py
	CXXFLAGS="$(CXXFLAGS) $(CXXWARNINGFLAGS)" CC="$(CC)" CXX="$(CXX)" $(PYTHON) -m build --wheel
	rm -f LICENSE.TXT
	touch .build-stamp

check-local: .build-stamp
	tox run

wheel: .build-stamp

clean-local:
	rm -rf build dist

endif
endif

src/QuantLib/quantlib_wrap.cpp src/QuantLib/QuantLib.py: ../SWIG/*.i
	$(SWIG) -python -c++ -outdir src/QuantLib -o src/QuantLib/quantlib_wrap.cpp ../SWIG/quantlib.i

dist-hook:
	mkdir -p $(distdir)/examples
	cp ./examples/*.py $(distdir)/examples
	mkdir -p $(distdir)/test
	cp ./test/*.py $(distdir)/test

EXTRA_DIST = README.txt pytest.ini setup.cfg setup.py swig.cmd tox.ini src/QuantLib/__init__.py $(BUILT_SOURCES)

