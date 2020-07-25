
Main changes for QuantLib-SWIG 1.19
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/12?closed=1>.

- We're sunsetting support for Python 2.7, which reached end of life
  in January 2020.  For the next couple of releases, we'll still check
  that the wrappers work with 2.7.  After that, we'll make no further
  effort to keep it working.

- Python examples can now be run as scripts as before, or as live
  notebooks on Binder.  They are available at
  <https://mybinder.org/v2/gh/lballabio/QuantLib-SWIG/binder?filepath=Python%2Fexamples>.

- Exported dividend barrier option and related engines.

- Exported choice of discretization for Heston process (thanks to
  GitHub user `feribg`).

- Added displacement parameter in `BlackCapFloorEngine` (thanks to
  Ralf Konrad).

- Exported Heston engine based on exponentially-fitted Laguerre
  quadrature rule (thanks to Klaus Spanderen).

- Exported spread options and Kirk spread option engine (thanks to
  Gorazd Brumen).

- Exported `AnalyticBSMHullWhiteEngine` class.

- Exported choice of timing adjustment for `BlackIborCouponPricer`
  (thanks to Matthias Lungwitz).

- Exported method for previous and next cash flow in `CashFlows`
  class.

- Detect correct location of include files when compiling C# wrappers
  via `make` (thanks to Ari Cooperman).

- Added support for VS 2019 in the solution for C# wrappers (thanks to
  Ralf Konrad).
