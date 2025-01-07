Main changes for QuantLib-SWIG 1.37
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/30?closed=1>.

- **Removed** the deprecated `SampledCurve` and `FixedRateBondForward`
  classes no longer available in the underlying C++ library;

- **Removed** the deprecated overload for `yoyInflationLeg`;

- Exported a number of new engines for basket and spread options;
  thanks to Klaus Spanderen (@klausspanderen).

- Exported Choi engine for Asian options; thanks to Klaus Spanderen
  (@klausspanderen).

- Exported new parameters and methods for `SwapRateHelper` and
  `OISRateHelper`; thanks to Eugene Toder (@eltoder) and Sotirios
  Papathanasopoulos (@sophistis42).

- Exported `MultipleResetsCoupon` and `MultipleResetsLeg` classes (@lballabio).

- Exported new constructors for `FittedBondDiscountCurve` (@lballabio).

- Exported additional arguments for `AssetSwap` constructor (@lballabio).

- Exported Wellington and Auckland variants for New Zealand calendar (@lballabio).

- Exported new constructors for YoY inflation curves (@lballabio).

- Exported KOFR index (@lballabio).

- Exported range-accrual coupon (@lballabio).

