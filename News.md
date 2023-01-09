
Main changes for QuantLib-SWIG 1.29
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/22?closed=1>.

- Enabled autodoc feature in Python; exported methods and classes have
  now docstrings reporting their interface and the types of the parameters.

- Enabled CI build and tests for the R wrappers; thanks to @AndLLA.

- **Removed** deprecated features no longer available in the
  underlying C++ library:
  - the constructor of `UnitedStates` missing an explicit market;
  - the `nominalTermStructure` method of `InflationTermStructure`;
  - the `CrossCurrencyBasisSwapRateHelper` class.

- Added `compounding` and `compoundingFrequency` parameters to
  `FixedRateLeg` (@lballabio).

- Exported `CashFlows::npvbps` method (@lballabio).

- Exported `baseFixing` and `indexFixing` methods in `IndexedCashFlow`
  (@lballabio).

- Exported new constructors for zero-inflation indexes (@lballabio).

- Exported missing arguments in `CreditDefaultSwap` constructor (@lballabio).

- Exported `Nearest` business-day convention (@lballabio).

- Exported `AmortizingCmsRateBond`; thanks to @chenyanlann.

- Exported `QuantoBarrierOption` and `QuantoBarrierEngine`; thanks to
  @chenyanlann.

- Avoided out-of-bound access to `Matrix` elements (@lballabio).

- Exported a number of LMM-related classes (@lballabio).

- Exported YoY inflation coupons and related classes (@lballabio).

- Exported the `CPI::laggedFixing` method; thanks to Marcin Rybacki
  (@marcin-rybacki).

- Exported `QdPlusAmericanEngine`, `QdFpAmericanEngine` and related
  classes; thanks to Klaus Spanderen (@klausspanderen).

- Added Python test case for Andreasen-Huge local volatility; thanks
  to Klaus Spanderen (@klausspanderen).

