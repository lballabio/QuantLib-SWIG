
CLEANFILES = cpp/quantlib_wrap.* cpp/*.so cpp/*.dylib csharp/*.cs *.dll

BUILT_SOURCES = cpp/quantlib_wrap.cpp cpp/quantlib_wrap.h

if HAVE_DOTNET
if BUILD_DOTNET

all-local: cpp/libNQuantLibc.@SHARED_LIB_EXTENSION@ csharp/bin/Release/net6.0/NQuantLib.dll

cpp/libNQuantLibc.@SHARED_LIB_EXTENSION@: cpp/quantlib_wrap.o
	$(CXX) -shared cpp/quantlib_wrap.o -o cpp/libNQuantLibc.@SHARED_LIB_EXTENSION@ `quantlib-config --libs`

cpp/quantlib_wrap.o: $(BUILT_SOURCES)
	$(CXX) -c -fpic $(CXXFLAGS) cpp/quantlib_wrap.cpp -o cpp/quantlib_wrap.o `quantlib-config --cflags`

csharp/bin/Release/net6.0/NQuantLib.dll: $(BUILT_SOURCES)
	$(DOTNET) build --nologo -c Release -p:Version=$(PACKAGE_VERSION) csharp/NQuantLib.csproj

nupkg: cpp/libNQuantLibc.@SHARED_LIB_EXTENSION@ csharp/bin/Release/net6.0/NQuantLib.dll
	$(DOTNET) pack --no-build -c Release -p:PackageVersion=$(PACKAGE_VERSION) --include-symbols --include-source csharp/NQuantLib.csproj

check-local: cpp/libNQuantLibc.@SHARED_LIB_EXTENSION@ csharp/bin/Release/net6.0/NQuantLib.dll
	ln -f cpp/libNQuantLibc.@SHARED_LIB_EXTENSION@ examples/
	cd examples && LD_LIBRARY_PATH=. $(DOTNET) run -c Release --project BermudanSwaption/BermudanSwaption.csproj
	cd examples && LD_LIBRARY_PATH=. $(DOTNET) run -c Release --project EquityOption/EquityOption.csproj
	cd examples && LD_LIBRARY_PATH=. $(DOTNET) run -c Release --project FiniteDifferenceMethods/FiniteDifferenceMethods.csproj
	cd examples && LD_LIBRARY_PATH=. $(DOTNET) run -c Release --project Times/Times.csproj
	rm -f examples/libNQuantLibc.@SHARED_LIB_EXTENSION@

clean-local:
	rm -rf csharp/bin csharp/obj examples/*/bin examples/*/obj

endif
endif

$(BUILT_SOURCES): ../SWIG/*.i
	$(SWIG) -csharp -c++ -outdir csharp \
            -namespace QuantLib -o cpp/quantlib_wrap.cpp ../SWIG/quantlib.i

dist-hook: $(BUILT_SOURCES)
	mkdir -p $(distdir)/cpp
	cp ./cpp/*.vcxproj $(distdir)/cpp
	cp ./cpp/quantlib_wrap.cpp $(distdir)/cpp
	cp ./cpp/quantlib_wrap.h $(distdir)/cpp
	mkdir -p $(distdir)/csharp
	cp ./csharp/*.csproj $(distdir)/csharp
	cp ./csharp/*.cs $(distdir)/csharp
	mkdir -p $(distdir)/examples
	mkdir -p $(distdir)/examples/BermudanSwaption
	cp ./examples/BermudanSwaption/BermudanSwaption.csproj $(distdir)/examples/BermudanSwaption
	cp ./examples/BermudanSwaption/BermudanSwaption.cs $(distdir)/examples/BermudanSwaption
	mkdir -p $(distdir)/examples/EquityOption
	cp ./examples/EquityOption/EquityOption.csproj $(distdir)/examples/EquityOption
	cp ./examples/EquityOption/EquityOption.cs $(distdir)/examples/EquityOption
	mkdir -p $(distdir)/examples/FiniteDifferenceMethods
	cp ./examples/FiniteDifferenceMethods/FiniteDifferenceMethods.csproj $(distdir)/examples/FiniteDifferenceMethods
	cp ./examples/FiniteDifferenceMethods/FiniteDifferenceMethods.cs $(distdir)/examples/FiniteDifferenceMethods
	mkdir -p $(distdir)/examples/Times
	cp ./examples/Times/Times.csproj $(distdir)/examples/Times
	cp ./examples/Times/Times.cs $(distdir)/examples/Times

EXTRA_DIST = \
    QuantLib.sln \
    QuantLib.props \
    swig.cmd \
    README.txt

