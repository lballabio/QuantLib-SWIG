
Main changes for QuantLib-SWIG 1.23
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/16?closed=1>.

- Exported overloaded constructors for piecewise inflation curves.

- Exported new `ZeroInflationCashFlow` class.

- Exported new constructor for `Currency` class (thanks to Marcin Rybacki).

- Exported `ZeroCouponSwap` class (thanks to Marcin Rybacki).

- Exported `MCDigitalEngine` class (thanks to Jack Gillett).

- Export updated 30/360 enumeration and constructors.

- Export `AnalyticHestonHullWhiteEngine`, `AnalyticH1HWEngine` and
  `FdHestonHullWhiteVanillaEngine` classes (thanks to Klaus Spanderen).

- Added payment lag and payment constructor to a few leg constructors
  (thanks to Marcin Rybacki).

- The `Type` enumeration defined in several swap classes was moved to
  their base `Swap` class.

- Updated ISDA CDS example in Python.  The differences between its
  results and Markit values are now within the desired tolerance
  (thanks to Francis Duffy).

- Removed constructors of piecewise yield and default curves taking an
  accuracy parameter (they were removed from the C++ library).

- Bond helper constructors now take a `BondPrice::Type priceType`
  argument instead of a `bool useCleanPrice` (the latter was removed
  from the C++ library).
