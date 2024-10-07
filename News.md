Main changes for QuantLib-SWIG 1.36
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/29?closed=1>.

- We're now using modern tooling to build and test Python wheels.
  Building now requires build besides setuptools, and testing requires
  pytest and tox. All of them can be installed in a virtual environment.

- **Removed** the deprecated constructors of the `ForwardRateAgreement` class.

- **Removed** the deprecated constructor of `YoYInflationIndex` taking
  a `ratio` parameter.

- **Removed** the deprecated `YYEUHICPr`, `YYFRHICPr`, `YYUKRPIr`,
  `YYUSCPIr` and `YYZACPIr` indexes.

- **Removed** the deprecated constructors of `CPICoupon` taking a
  `spread` parameter and its `spread` method, as well as the
  deprecated `withSpreads` method of `CPILeg`.

- **Breaking**: in Python, the multiplication between two `ql.Array`
  instances would return the dot product.  It now returns the
  element-wise product, like in C++.  Also, exposed more operators.
  Thanks to Eugene Toder (@eltoder).

- Exported `SpreadedSwaptionVolatility` class (@lballabio).

- Exported `Index::pastFixing` and the constructor of `EquityIndex`
  taking currency information; thanks to Ralf Konrad Eckel
  (@ralfkonrad).

- Exported specialized Warsaw Stock Exchange (WSE) calendar for
  Poland; thanks to Marcin Bogusz (@marcinfair).

- Exported missing volatility-type parameter for SABR interpolation
  (@lballabio).  This allows using it for normal volatilities.

- Exported `startOfMonth` and `isStartOfMonth` methods for both `Date`
  and `Calendar` (@lballabio).

- Exported `CompoundingOvernightIndexedCouponPricer` and
  `ArithmeticAveragedOvernightIndexedCouponPricer`, and export
  corresponding pricer parameter for the `OISRateHelper` and
  `DatedOISRateHelper` constructors (@lballabio).

- Export additional custom-constraint parameter for non-linear fitting
  methods (@lballabio).

- Exported `needsForecast` and `lastFixingDate` methods for inflation
  indexes (@lballabio).

- Exported new optimizer and end-criteria parameters for the
  `GlobalBootstrap` constructor (@lballabio).

- Exported new interpolation parameter for YoY inflation coupons
  (@lballabio).
