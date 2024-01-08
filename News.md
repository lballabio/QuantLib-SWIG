
Main changes for QuantLib-SWIG 1.33
===================================

More details on the changes are available in ChangeLog.txt and at
<https://github.com/lballabio/QuantLib-SWIG/milestone/26?closed=1>.

- Exported Burley 2020 Sobol generator (@lballabio).

- Allowed different calendars and frequencies for different legs in
  `OISRateHelper`; thanks to Eugene Toder (@eltoder).

- Exported convex-monotone forward-rate curve (@lballabio).

- Exported support for angled contour shift integrals in Heston model;
  thanks to Klaus Spanderen (@klausspanderen).

- Allowed negative payment lag in swap legs; thanks to Fredrik Gerdin
  Börjesson (@gbfredrik).

- Exported `reset` method in calendars; thanks to Fredrik Gerdin
  Börjesson (@gbfredrik).

- Added Python tests for `BondFunctions`; thanks to Francois Botha
  (@igitur).
