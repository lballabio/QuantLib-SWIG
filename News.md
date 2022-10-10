
Main changes for QuantLib-SWIG 1.28
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/21?closed=1>.

- **Removed** deprecated features no longer available in the
  underlying C++ library:
  - the constructors of `ZeroCouponInflationSwap` and
    `ZeroCouponInflationSwapHelper` missing an explicit CPI
    interpolation type;
  - the constructors of `ActualActual` and `Thirty360` missing an
    explicit choice of convention.

- **Renamed** `RelinkableYoYOptionletVolatilitySurface` to
  `RelinkableYoYOptionletVolatilitySurfaceHandle`.  The old name is
  still available in Python as deprecated. Currently we have no way to
  do so in other languages.

- Added an implicit conversion in C# from `bool` to
  `boost::optional<bool>`, making it possible to pass parameters of
  this type.  Python already had typemaps defined.  Other languages
  can pass `OptionalBool(b)` where `b` is the desired bool.

- Exported the `Gaussian1dCapFloorEngine` class; thanks to @jacek-bator.

- Exported `LazyObject` methods in `PiecewiseYieldCurve`; thanks to
  Francois Botha (@igitur).

- Exported Act/366 and Act/365.25 day counters; thanks to Ignacio
  Anguita (@IgnacioAnguita).

- Exported `PartialTimeBarrierOption` class and related engine; thanks
  to Ignacio Anguita (@IgnacioAnguita).

- Added missing `operator-` to `Date` in C#.

- Added a few default parameters to the `SABRInterpolation` constructor.

- Exported new constructor for `SabrSmileSection`.

- Exported new `sinkingSchedule` and `sinkingNotionals` functions.

- Exported new overload for `CallableBond::impliedVolatility`.

- Exported missing end-of-month optional parameter for `OISRateHelper`
  constructor.

