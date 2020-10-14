
Main changes for QuantLib-SWIG 1.20
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/13?closed=1>.

- We're sunsetting support for Python 2.7, which reached end of life
  in January 2020.  For the next release, we'll still check that the
  wrappers work with 2.7.  After the next release, we'll make no
  further effort to keep it working.

- SWIG wrappers now work also if the C++ library was compiled using
  `std::shared_ptr` instead of `boost::shared_ptr` (thanks to Joseph
  Wang).

- The `BaroneAdesiWhaleyApproximationEngine` and
  `BjerksundStenslandApproximationEngine` classes used to be renamed
  to `BaroneAdesiWhaleyEngine` and `BjerksundStenslandEngine`,
  respectively.  This is no longer the case.

- Exported mixing factor to Heston SLV process and engines (thanks to
  Jack Gillett).

- Exported a number of inflation-related classes (thanks to Matthias
  Lungwitz).

- Exported Crank-Nicolson finite-differences scheme (thanks to Klaus
  Spanderen).

- Exported `SwaptionVolatilityCube` class (thanks to Marcin Rybacki).

- Exported Cox-Ingersoll-Ross short-rate model.

- Exported callable zero-coupon bond.

- Exported SABR interpolation.

- Made the `Date` class comparable and convertible to string in C#.

