Changes for QuantLib-SWIG 1.41
==============================

Removed features
----------------

Features removed from the C++ library in this release were also
removed from these wrappers; see
<https://github.com/lballabio/QuantLib-SWIG/pull/783> for a full list.

Possibly breaking changes
-------------------------

- Due to limitations of SWIG in case of overloaded methods, the
  constructors of `ContinuousAveragingAsianOption` and
  `ContinuousArithmeticAsianLevyEngine` no longer accept keyword
  arguments in Python;
- In the constructor of `ZeroCouponInflationSwapHelper`, the
  incorrectly-named `bcd` argument was renamed to `bdc` (for
  business-day convention); this can break calls using that keyword
  argument.


Full list of pull requests
--------------------------

All the pull requests merged in this release are listed on its release
page at <https://github.com/lballabio/QuantLib-SWIG/releases/tag/v1.41>.

The list of commits since the previous release is available in `ChangeLog.txt`.
