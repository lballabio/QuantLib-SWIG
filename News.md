Changes for QuantLib-SWIG 1.42
==============================

Starting from this release, the Python wrappers provide better support for the new free-threading interpreter in Python 3.14 (the ones in version 1.41 worked but would not allow the interpreter to disable the GIL).

Removed features
----------------

Features removed from the C++ library in this release were also removed from these wrappers; see <https://github.com/lballabio/QuantLib-SWIG/pull/813> for a full list.


Full list of pull requests
--------------------------

All the pull requests merged in this release are listed on its release page at <https://github.com/lballabio/QuantLib-SWIG/releases/tag/v1.42>.

The list of commits since the previous release is available in `ChangeLog.txt`.
