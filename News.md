
Main changes for QuantLib-SWIG 1.24
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/17?closed=1>.

- Breaking change: removed inflation-curve constructors taking a
  nominal curve (they were deprecated and were removed from the C++
  library in version 1.24).

- Breaking change: removed the long-deprecated
  `BaroneAdesiWhaleyEngine` and `BjerksundStenslandEngine` aliases for
  `BaroneAdesiWhaleyApproximationEngine` and
  `BjerksundStenslandApproximationEngine`, respectively.

- Exported `CliquetOption` class and related pricing engines (thanks
  to Jack Gillett).

- Made the `Period` class equatable and comparable in C# (thanks to
  Ralf Konrad).

- Exported the missing `endOfMonth` parameter in `SwapRateHelper`
  constructor (thanks to Fidel Selva).

- Exported the new `ConstNotionalCrossCurrencyBasisSwapRateHelper` and
  `MtMCrossCurrencyBasisSwapRateHelper` rate helpers (thanks to Marcin
  Rybacki).

- Exported the new `RiskyBondEngine` class.

- Exported the new Chilean calendar.

- Exported the new `ThirdWednesdayInclusive` date-generation rule.

- Exported the new `useIndexedCoupon` parameter in the constructors of
  `BlackIborCouponPricer`, `IborLeg`, `SwapRateHelper` and
  `VanillaSwap`.
