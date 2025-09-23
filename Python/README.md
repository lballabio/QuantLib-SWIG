
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
your QuantLib directory (e.g., `C:\Lib\QuantLib`.)

Once built, the resulting wheel can be installed with `pip`.

Finally, testing the wheel requires `tox`, also available via
`pip install`.  Once available, running `tox run` will run the
test suite.

