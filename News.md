
Main changes for QuantLib-SWIG 1.16
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/9?closed=1>.

Starting from this release, the SWIG wrappers make use of the native
SWIG support for `shared_ptr`.  This caused us to drop support for
Perl, which doesn't provide the functionality.  If you're interested
in keeping the Perl module alive, consider contacting the SWIG
maintainers and asking advice about providing native `shared_ptr`
support for that language.

New features exported in this release:
- CEV engines (thanks to Klaus Spanderen).
- more methods of the Schedule class (thanks to Matthias Lungwitz).
- Ju quadratic engine (thanks to Matthias Lungwitz).
- more Ibor indexes (thanks to Matthias Lungwitz).
- more results for caps and floors (thanks to Wojciech Slusarski).
- CDS options and Black engine (thanks to Matthias Lungwitz).
- Andreasen huge volatility (thanks to Matthias Lungwitz).
- Overnight-indexed swap index (thanks to Matthias Lungwitz).
- Kirk engine for basket options (thanks to Matthias Lungwitz).
- G2 processes.
