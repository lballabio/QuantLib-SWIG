
Main changes for QuantLib-SWIG 1.31
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/24?closed=1>.

- **Removed** deprecated features no longer available in the
  underlying C++ library:
  - The `CPICoupon` constructor taking a number of fixing days and its
    `adjustedFixing` method.
  - The `withFixingDays` methods of `CPILeg`.
  - The `ZeroInflationCashFlow` constructor taking a calendar and
    business-day convention.
  - The `LexicographicalView` class.

- Exported new U.S. SOFR calendar (@lballabio).

- Exported new constructors and `indexRatio` method for `CPICoupon`
  (@lballabio).

- Exported new constructors and `underlyingIndex` method for
  `YoYInflationIndex` (@lballabio).

- Exported new constructors for `ForwardRateAgreement` (@lballabio).

- Rework Python tests to follow standard conventions; thanks to Eugene
  Toder (@eltoder).

- Updated constructor of `DatedOISRateHelper` to take new parameters;
  thanks to Eugene Toder (@eltoder).

- Exported missing currencies and crypto; thanks to Fredrik Gerdin
  Börjesson (@gbfredrik).

- Exported `LogMixedLinearCubic` interpolator and corresponding
  discount curves; thanks to Eugene Toder (@eltoder).

- Exported `ArithmeticAverageOIS` and the corresponding rate helper;
  thanks to Eugene Toder (@eltoder).

- Exported a few missing inspectors for `Swap`; thanks to Eugene Toder
  (@eltoder).

- Exported CORRA, SWESTR and DESTR indexes; thanks to Fredrik Gerdin
  Börjesson (@gbfredrik).

- Exported new constructor and Python tests for `JointCalendar`;
  thanks to Fredrik Gerdin Börjesson (@gbfredrik).

- Exported new LazyObject interface (@lballabio).

- Added Python examples for callable bonds and caps; thanks to Nijaz
  Kovacevic (@NijazK).

- Added convenience methods `of` and `toLocalDate` to Java wrappers
  that convert QuantLib dates from and to `java.time.LocalDate`; and
  example is provided.  Thanks to Ralf Konrad (@ralfkonrad).
