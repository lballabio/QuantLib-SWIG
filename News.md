
Main changes for QuantLib-SWIG 1.17
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/10?closed=1>.

- The Python module is now available on PyPI as simply `QuantLib`.
  The old name `QuantLib-Python` will still be provided for a while
  but it will cause the new module to be installed.

- As of this release, the Ruby wrappers are deprecated.  They will be
  removed in next release.  They have been broken for a while, and
  nobody expressed any interest in fixing them.

- It is now possible to specify at run time whether to use indexed
  coupons (thanks to Ralf Konrad).  This makes it possible to enable
  them using prebuilt binaries without having to recompile them.  You
  can call either of the static methods `IborCoupon::createAtParCoupons` or
  `IborCoupon::createIndexedCoupons` to specify your preference.  For
  the time being, the methods above must necessarily be called before
  creating any instance of `IborCoupon` or of its derived classes.

- Exported `MakeOIS` class (thanks to Wojciech Slusarski).

- Exported American quanto dividend option engine (thanks to Klaus
  Spanderen).

- Exported `legBPS` method for swaps (thanks to Mike DelMedico).

- Exported cubic spline variant of piecewise discount curve,
  interpolated discount curve and interpolated zero-rate curve (thanks
  to Mike DelMedico).

- Exported interpolated zero and YoY inflation curves (thanks to
  Matthias Lungwitz).

- Exported overnight-indexed coupon and leg (thanks to Matthias
  Lungwitz).

- Exported a few overnight LIBOR indexes.

- Exposed ex-coupon functionality for floaters.

- Added static convenience methods to call different fixed-rate bond
  constructors using Python keyword arguments (thanks to Prasad
  Somwanshi).

- Exported Newton solver to Java and C++ (thanks to Klaus Spanderen).
