
Main changes for QuantLib-SWIG 1.32
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/25?closed=1>.

- Avoid using the deprecated `distutils` module for the Python
  wrappers; `setuptools` is now required for building (@lballabio).

- Exported `LastFixingQuote`; thanks to Eugene Toder (@eltoder).

- Added `redemptions` and `paymentLag` arguments to amortizing bond
  constructors; thanks to Gyan Sinha (@gyansinha).

- Exported utility function to simplify notification graph (@lballabio).

- Exported a few exotic options (Margrabe, compound, chooser) and
  related engines (@lballabio).

- Exported new constructor for OIS (@lballabio).

- Exported missing parameters for iterative bootstrap (@lballabio).

- Exported Xoshiro256** RNG (@lballabio).

