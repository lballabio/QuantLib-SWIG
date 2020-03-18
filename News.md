
Main changes for QuantLib-SWIG 1.18
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/11?closed=1>.

- As announced in the past release, the Ruby wrappers were removed.
  They have been broken for a while, and nobody expressed any interest
  in fixing them.

- Exported most of the inner machinery (meshers, operators, boundary
  conditions, schemes, solvers...) of the finite-difference framework
  (thanks to Klaus Spanderen).

- Exported GJR-GARCH process, model, analytic engine and MC engine
  (thanks to Pedro Coelho).

- The accuracy of piecewise curve can now be passed as an argument
  to the `IterativeBootstrap` class, which in turn can be passed to
  the curve.  The new class also allows to set minimum and maximum
  values explicitly.

- Exported the new `GlobalBootstrap` class and the corresponding
  `GlobalLinearSimpleZeroCurve` curve.

- Exported the `CmsMarket` class (thanks to Matthias Lungwitz).

- Exported convex monotone and Kruger cubic and log-cubic
  interpolation (thanks to Miguel Villasmil).

- Exported `InflationCoupon` and `CPICoupon` classes with
  corresponding functions `as_inflation_coupon` and `as_cpi_coupon`.

- Exported missing methods of the `SwaptionVolatilityStructure` class
  (thanks to Matthias Lungwitz).

- Exported the `CallableFixedRateBond` class and a few missing methods
  of the `CallableBond` class.

- Exported the `enforcesTodaysHistoricFixings` flag from the
  `Settings` class (thanks to Tomáš Křehlík).

- Exported the `OvernightIndexFutureRateHelper` class (thanks to
  Miguel Villasmil).

- Exported the `SofrFutureRateHelper` class.

- Allowed use of normal volatility with the `CapHelper` class.

