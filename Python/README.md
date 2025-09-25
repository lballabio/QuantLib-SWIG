
The C++ wrappers for the QuantLib-Python extension module are created
by means of SWIG (Simplified Wrapper and Interface Generator)
available from <http://www.swig.org/>; the latest version is suggested.

Building the wrappers requires the `setuptools` and `build` packages;
both can be installed (preferably in a virtual environment) by means
of `pip install`.

Generating the wrappers is not required if you are using a distributed
tarball. If you're building from a Git checkout, instead, use the
command `swig.cmd` on Windows or `make` on other platforms.  Running
`make` also builds the wrappers as a wheel; on Windows, this requires
the explicit command `python -m build --wheel` instead.

The build step requires that the QuantLib headers and library can be
found by the compiler. On Unix-like platforms, this requires that
`quantlib-config` is in your path. On the Windows platform, instead,
it requires you to define a `QL_DIR` environment variable pointing to
your QuantLib directory (e.g., `C:\Lib\QuantLib`.) Another environment
variable `QL_DEBUG` on Windows should be set to `TRUE` if you are building
against a debug version of QuantLib with the objective to
debug the library called from Python. On Unit-like platforms, appropriate flags
can be setup in `CXXFLAGS` and `LDFLAGS`. If boost is required for compilation,
it can be set via the environment variable `INCLUDE`.


Once built, the resulting wheel can be installed with `pip`.

Finally, testing the wheel requires `tox`, also available via
`pip install`.  Once available, running `tox run` will run the
test suite.

