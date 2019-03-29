
Main changes for QuantLib-SWIG 1.15
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/8?closed=1>.

- Add `settlementMethod` parameter to `Swaption` constructor.

- Add `polynomOrder`, `polynomType` and `nCalibrationSamples`
  parameters to `MCAmericanBasketEngine` constructor.


Starting from next release, the SWIG wrappers will make use of the
native SWIG support for `shared_ptr`.  This will cause us to drop
support for Perl, which doesn't provide the functionality.  If you're
interested in keeping the Perl module alive, consider contacting the
SWIG maintainers and asking advice about providing native `shared_ptr`
support for that language.
