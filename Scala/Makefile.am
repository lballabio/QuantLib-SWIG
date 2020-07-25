
EXAMPLES = CPIBond EquityOptions HestonModelCalibration HestonMonteCarlo ObserverPattern Optimizer RandomNumbers Swing

SWIGFLAGS =

if JAVA_AUTOLOAD
SWIGFLAGS += -DJAVA_AUTOLOAD
endif

if JAVA_FINALIZER
SWIGFLAGS += -DJAVA_FINALIZER
endif

if JAVA_AUTOCLOSEABLE
SWIGFLAGS += -DJAVA_AUTOCLOSEABLE
endif

if HAVE_JAVA
if BUILD_JAVA
if HAVE_SCALA
if BUILD_SCALA

examples/%.class: examples/%.scala
	$(SCALAC) -cp ../Java/QuantLib.jar -language:postfixOps examples/$*.scala

.PHONY: $(EXAMPLES)

check-local: $(EXAMPLES)

$(EXAMPLES): $(EXAMPLES:%=examples/%.class)
	LD_LIBRARY_PATH=../Java $(SCALA) -cp ".:../Java/QuantLib.jar" examples.$@

EXTRA_DIST = README.txt

endif
endif
endif
endif


dist-hook:
	mkdir -p $(distdir)/examples
	cp ./examples/*.scala $(distdir)/examples

clean-local:
	rm -f examples/*.class
