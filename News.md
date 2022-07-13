
Main changes for QuantLib-SWIG 1.27
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/20?closed=1>.


- Fixed code generation when Java's `AutoCloseable` support is enabled
  through the `--enable-java-autocloseable` configure switch.

- Exported the `SviSmileSection` and `SviInterpolatedSmileSection`
  classes; thanks to Fredrik Gerdin Börjesson (@gbfredrik).

- Exported the `QuantoTermStructure` class; thanks to Sebastian Bohlen
  (@BohlSeb).

- Exported the `TurnbullWakemanAsianEngine` class; thanks to Jack
  Gillett (@jackgillett101).

- Exported shorter name (the same as in C++) for
  `PiecewiseZeroSpreadedTermStructure`.  The older and uglier
  `SpreadedLinearZeroInterpolatedTermStructure` is still available.

- Exported the `previousCashFlowAmount` and `nextCashFlowAmount`
  methods from the `CashFlows` class; thanks to Marcin Rybacki
  (@marcin-rybacki).

- Exported a few missing methods from the `CreditDefaultSwap` class.

- Removed the `FDShoutEngine` class, no longer available in the
  underlying C++ library; thanks to Fredrik Gerdin Börjesson
  (@gbfredrik).

- Removed reference to the deprecated `Disposable` class from
  interface files; thanks to Jonathan Sweemer (@sweemer).

