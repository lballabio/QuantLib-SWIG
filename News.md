
Main changes for QuantLib-SWIG 1.22
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/15?closed=1>.

- As previously announced, this release drops support for
  Python 2.7, which reached end of life in January 2020.

- Exported revised `SubPeriodCoupon` class (thanks to Marcin Rybacki).

- Exported new `FdBlackScholesShoutEngine` class.

- Exported optional discount curve in AnalyticEuropeanEngine constructor.

- Exported the `CrossCurrencyBasisSwapRateHelper` (thanks to Marcin Rybacki).

- Exported new constructors for Asian options (thanks to Jack Gillett).

- Exported new method `hasHistoricalFixing` for indexes (thanks to Ralf Konrad).

- Exported revised `OvernightIndexFuture` interface.

- Classes `CallabilityPrice`, `FDBermudanEngine`, `FDEuropeanEngine`,
  `FDAmericanEngine`, `FDDividendEuropeanEngine` and
  `FDDividendAmericanEngine` were removed in the C++ library after a
  deprecation period.

