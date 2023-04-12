
Main changes for QuantLib-SWIG 1.30
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/23?closed=1>.

- **Removed** deprecated features no longer available in the
  underlying C++ library:
  - the `WulinYongDoubleBarrierEngine` alias for `SuoWangDoubleBarrierEngine`;
  - the `spotIncome` and `spotValue` methods of `ForwardRateAgreement`;
  - constructors for `InterpolatedZeroInflationCurve` and
    `PiecewiseZeroInflationCurve` taking an `indexIsInterpolated`
    parameter;
  - the `indexIsInterpolated` method of `InflationTermStructure`;
  - some overloaded constructors of `SofrFutureRateHelper`.

- **Renamed** `SwaptionVolCube1` to `SabrSwaptionVolatilityCube` and
  `SwaptionVolCube2` to `InterpolatedSwaptionVolatilityCube`, as in
  the underlying C++ library; the old names remain available in Python
  but not in other languages.

- Exported new `EquityCashFlow`, `EquityIndex` and
  `EquityTotalReturnSwap` classes with a few tests; thanks to Marcin
  Rybacki (@marcin-rybacki).

- Exported constructors for vanilla and barrier pricing engines taking
  discrete dividends; this makes `DividendVanillaOption` and
  `DividendBarrierOption` obsolete (@lballabio).

- Exported new calendars for Austria, Botswana and Romania; thanks to
  Fredrik Gerdin Börjesson (@gbfredrik).

- Exported new ASX calendar for Australia (@lballabio).

- Exported `FixedLocalVolSurface` and `GridModelLocalVolSurface`
  classes with a test; thanks to Klaus Spanderen (@klausspanderen).

- Exported new CPICoupon constructors (@lballabio).

- Exported UKHICP index (@lballabio).

- Exported a few African currencies (@lballabio).

