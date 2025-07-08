Main changes for QuantLib-SWIG 1.39
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/32?closed=1>.

- **Removed** features deprecated in version 1.34 and no longer available
  in the underlying C++ library:
  - the overloads of `Bond::yield`, `BondFunctions::atmRate`,
    `BondFunctions::yield` and `BondFunctions::zSpread` taking a price
    as a `Real` instead of a `Bond::Price` instance;
  - the `Swaption::underlyingSwap` and
    `SwaptionHelper::underlyingSwap` methods;
  - the constructors of `InterpolatedZeroInflationCurve`,
    `InterpolatedYoYInflationCurve`, `PiecewiseZeroInflationCurve` and
    `PiecewiseYoYInflationCurve` taking an observation lag;
  - the overload of `InflationTermStructure::setSeasonality` taking no arguments;
  - the `fixedRateBond` method of the `FixedRateBondHelper` class.

- Added preliminary support for the new free-threading Python
  interpreter; thanks to Klaus Spanderen (@Klausspanderen).  No wheels
  are provided for it at this time.

- Java compilation flags can now be passed by setting the
  `JAVAC_FLAGS` environment variable; thanks to @UnitedMarsupials.

- Exported `convexityAdjustment` method for `FuturesRateHelper` and
  `OvernightIndexFutureRateHelper` classes; thanks to Eugene Toder
  (@eltoder).

- Passing a nominal curve to the `ZeroCouponInflationSwapHelper`
  constructor is now optional (@lballabio).

- The `OISRateHelper` constructor can now take a calendar for the
  overnight leg; thanks to Eugene Toder (@eltoder).

- Exported the `CustomIborIndex` class; thanks to Eugene Toder
  (@eltoder).

- Exported the `sabrGuess` function (@lballabio).

- Exported the `SARON` index (@lballabio).

- Exported the static `FxSwapRateHelper.forDates` method; thanks to
  Eugene Toder (@eltoder).

- The `OptionletStripper1` constructor can be passed a frequency so
  that it can be used with overnight indexes (@lballabio).

- Exported the SHIR calendar (@lballabio).
