
Main changes for QuantLib-SWIG 1.26
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/19?closed=1>.

- Running `make` in the `CSharp` folder (after `../configure` in the
  root folder) now uses `dotnet` instead of older C# compilers.
  The existing VC# projects were also updated to .Net.

- Fixed a compilation error in wrappers when using strict C++17 mode
  (thanks to Jonghee Lee).

- Exported a few more constructors for `FraRateHelper` (thanks to Marcin Rybacki).

- Exported the `SpreadFittingMethod` class.

- Exported the new `BondForward` class (thanks to Marcin Rybacki).

- Ensure that periods that compare equal have the same hash in C# and Python.

- Exported `SplineLogCubic` interpolation and the corresponding
  `NaturalLogCubicDiscountCurve`, `PiecewiseNaturalCubicZero` and
  `PiecewiseNaturalLogCubicDiscount` classes. Also exported
  `KrugerLogDiscountCurve` and `KrugerZeroCurve` based on Kruger
  interpolation (thanks to Marcin Rybacki).

- Exported `as_overnight_indexed_coupon` function to downcast such
  coupons when needed (thanks to Marcin Rybacki).

