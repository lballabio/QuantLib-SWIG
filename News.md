Main changes for QuantLib-SWIG 1.34
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/27?closed=1>.

- Upgrade to SWIG 4.2.x.  This allows to use Python's limited API and
  thus reduce the number of official wheels to cover the same Python
  versions.

- Allow swaptions to use OIS as underlying (@lballabio).

- Pass explicit base date to inflation curves instead of observation
  lag (@lballabio).

- Exported `SavedSettings` as a context manager in Python; thanks
  to Eugene Toder (@eltoder).

- Exported parabolic (Hermite) cubic spline interpolation schemes;
  thanks to Marcin Rybacki (@marcin-rybacki).

- Exported additional interpolation schemes for
  `InterpolatedPiecewiseZeroSpreadedTermStructure`; thanks to Marcin
  Rybacki (@marcin-rybacki).

- Exported Tona index; thanks to Jonghee Lee (@nistick21).

- Removed inflation index constructors with `interpolated` parameters
  as well as the `interpolated` method in `InflationIndex`.  They're
  no longer available in C++ (@lballabio).

- Export a few new methods for MakeOIS and MakeVanillaSwap; thanks to
  Eugene Toder (@eltoder).

- Exported `cdsMaturity` function (@lballabio).

- Enable different definition of macro `QL_JAVA_INTERFACES`; thanks to
  Ralf Konrad (@ralfkonrad).

- Define a few additional operators in C++ instead of Python; thanks
  to Eugene Toder (@eltoder).

- Removed uncallable internal `EndCriteria::operator()` method
  (@lballabio).
