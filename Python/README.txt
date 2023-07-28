
The C++ wrappers for the QuantLib-Python extension module are created
by means of SWIG (Simplified Wrapper and Interface Generator)
available from <http://www.swig.org/>; the latest version is suggested.

Building the wrappers requires the Setuptools package.

Generating the wrappers is not required if you are using a distributed
tarball. If you're building from a Git checkout, instead, use the
command:

    python setup.py wrap

The commands to be issued for building and testing the
wrappers are:

    python setup.py build
    python setup.py test

respectively.

The build step requires that the QuantLib headers and library can be
found by the compiler. On Unix-like platforms, this requires that
`quantlib-config` is in your path. On the Windows platform, instead,
it requires you to define a `QL_DIR` environment variable pointing to
your QuantLib directory (e.g., `C:\Lib\QuantLib`.)

The suggested way to install the module is to run first

    python setup.py bdist_wheel

which requires the wheel module to be available and which will create
a QuantLib wheel in the dist folder.  The resulting wheel can then
be installed with pip in the desired environment.

The test suite is implemented on top of the PyUnit framework, which is
included in the Python standard library.
