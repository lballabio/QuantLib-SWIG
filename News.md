
Main changes for QuantLib-SWIG 1.21
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/14?closed=1>.

- As previously announced, this release is the last one to support
  Python 2.7, which reached end of life in January 2020.

- Exported methods for `VolDeltaQuote` class (thanks to Jack Gillett).

- Exported `localVolatility()` method for Black-Scholes process
  (thanks to Jack Gillett).

- Exported `type()` method for vanilla swaps (thanks to Ralf Konrad).

- Exported constructors with full parameter lists for CDS helpers
  (thanks to Joe Song).

- Exported amortizing-bond constructor taking an `InterestRate`
  instance (thanks to Piter Dias).

- Exported Sobol-based multi-path generator (thanks to Jack Gillett).

- Exported Monte Carlo and analytic forward option engines based on
  the Heston model (thanks to Jack Gillett).

- Exported ultimate-forward term structure (thanks to Marcin Rybacki).

- Exported a few Monte Carlo and analytic Asian option engines based
  on the Heston model (thanks to Jack Gillett).

- Exported swap constructor taking multiple legs.

- Exported lookback options.

- Exported overnight-index futures.

- Avoided memory access issue with path generators in Python (thanks
  to Klaus Spanderen for the heads-up).

- Added an example of cash-flow analysis in Python.
