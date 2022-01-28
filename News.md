
Main changes for QuantLib-SWIG 1.25
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/18?closed=1>.

- **Breaking change:** exported updated interface for convertible bonds and their engine.

- **Breaking change (except for Python):** renamed `WulinYongDoubleBarrierEngine`
  to `SuoWangDoubleBarrierEngine`.

- Added a few missing methods to `Schedule` (thanks to Ralf Konrad).

- Exported `CPICoupon`, `CPICashFlow`, `CPILeg`.

- Exported new argument to `SabrSmileSection` constructor to allow normal volatilities.

- Exported new constructor and `amount` method for `ForwardRateAgreement`.

- Exported new constructors for `SofrFutureRateHelper`.

- Exported new constructors for zero-inflation curves.

- Exported a few more finite-difference classes (thanks to Klaus Spanderen).

- Exported new basis-swap rate helpers.

- Exported `ESTR` class (thanks to Kirill Egorov).

- Exported `StrippedOptionlet` class.
