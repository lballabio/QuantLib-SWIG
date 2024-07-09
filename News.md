Main changes for QuantLib-SWIG 1.35
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/28?closed=1>.

- Removed deprecated classes `DividendVanillaOption` and
  `DividendBarrierOption`.
  
- Removed deprecated constructor of `AnalyticDividendEuropeanEngine`
  taking only a process and no dividends.

- Exported missing `CashAnnuityModel` parameter for Black and
  Bachelier swaption engines (@lballabio).

- Exported Ziggurat Gaussian RNG; thanks to Ralf Konrad Eckel
  (@ralfkonrad).

- Exported a few missing `CashFlows` methods (@lballabio); thanks to
  GitHub user @heiieh for the heads-up.

- Exported new `IborCoupon::hasFixed` method (@lballabio).

- Exported new `FittedBondDiscountCurve::resetGuess` method (@lballabio).

- `EuriborSW` renamed to `Euribor1W`, old name still available for a
  while (@lballabio).

- Exported lookback days, lockout days and observation shift for
  overnight-indexed coupons, swaps and helpers (@lballabio).

- Exported `SimpleQuote::reset` method; thanks to Eugene Toder
  (@eltoder).

