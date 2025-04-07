Main changes for QuantLib-SWIG 1.38
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/31?closed=1>.

- **Removed** the deprecated `Currency` constructor no longer
  available in the underlying C++ library;

- Exported forward curve with a number of additional interpolations;
  thanks to Sotirios Papathanasopoulos (@sophistis42) and to
  @paolodelia99.

- Exported `FuturesConvAdjustmentQuote`; thanks to Eugene Toder
  (@eltoder).

- Exported missing default parameters for `MakeVanillaSwap` and
  `MakeOIS`; thanks to Eugene Toder (@eltoder).

- Exported new constructors for `DepositRateHelper` and
  `FraRateHelper`; thanks to Eugene Toder (@eltoder).

- Exported new constructor arguments for cross-currency basis-swap
  helpers; thanks to @kp9991-git.

- Exported methods to return the underlying process from a few models
  (@lballabio).

- Exported new constructors for YoY inflation indexes (@lballabio).

- Exported a few more exotic options and engines (@lballabio):
  - `TwoAssetBarrierOption` with `AnalyticTwoAssetBarrierEngine`;
  - `HolderExtensibleOption` with `AnalyticHolderExtensibleOptionEngine`;
  - `WriterExtensibleOption` with `AnalyticWriterExtensibleOptionEngine`;
  - `TwoAssetCorrelationOption` with `AnalyticTwoAssetCorrelationEngine`;
  - `AnalyticPDFHestonEngine`.

- Exported piecewise forward-spreaded term structure (@lballabio).
