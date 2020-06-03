---
jupyter:
  jupytext:
    formats: py:percent,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.4.2
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Interest-rate swaps

Copyright (&copy;) 2004, 2005, 2006, 2007 StatPro Italia srl

This file is part of QuantLib, a free-software/open-source library
for financial quantitative analysts and developers - https://www.quantlib.org/

QuantLib is free software: you can redistribute it and/or modify it under the
terms of the QuantLib license.  You should have received a copy of the
license along with this program; if not, please email
<quantlib-dev@lists.sf.net>. The license is also available online at
<https://www.quantlib.org/license.shtml>.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the license for more details.

```python
import QuantLib as ql
import pandas as pd
```

### Global data

```python
calendar = ql.TARGET()
todaysDate = ql.Date(6, ql.November, 2001)
ql.Settings.instance().evaluationDate = todaysDate
settlementDate = ql.Date(8, ql.November, 2001)
```

### Market quotes

```python
deposits = {
    (3, ql.Months): 0.0363,
}
```

```python
FRAs = {(3, 6): 0.037125, (6, 9): 0.037125, (9, 12): 0.037125}
```

```python
futures = {
    ql.Date(19, 12, 2001): 96.2875,
    ql.Date(20, 3, 2002): 96.7875,
    ql.Date(19, 6, 2002): 96.9875,
    ql.Date(18, 9, 2002): 96.6875,
    ql.Date(18, 12, 2002): 96.4875,
    ql.Date(19, 3, 2003): 96.3875,
    ql.Date(18, 6, 2003): 96.2875,
    ql.Date(17, 9, 2003): 96.0875,
}
```

```python
swaps = {
    (2, ql.Years): 0.037125,
    (3, ql.Years): 0.0398,
    (5, ql.Years): 0.0443,
    (10, ql.Years): 0.05165,
    (15, ql.Years): 0.055175,
}
```

We'll convert them to `Quote` objects...

```python
for n, unit in deposits.keys():
    deposits[(n, unit)] = ql.SimpleQuote(deposits[(n, unit)])
for n, m in FRAs.keys():
    FRAs[(n, m)] = ql.SimpleQuote(FRAs[(n, m)])
for d in futures.keys():
    futures[d] = ql.SimpleQuote(futures[d])
for n, unit in swaps.keys():
    swaps[(n, unit)] = ql.SimpleQuote(swaps[(n, unit)])
```

...and build rate helpers.

```python
dayCounter = ql.Actual360()
settlementDays = 2
depositHelpers = [
    ql.DepositRateHelper(
        ql.QuoteHandle(deposits[(n, unit)]),
        ql.Period(n, unit),
        settlementDays,
        calendar,
        ql.ModifiedFollowing,
        False,
        dayCounter,
    )
    for n, unit in deposits.keys()
]
```

```python
dayCounter = ql.Actual360()
settlementDays = 2
fraHelpers = [
    ql.FraRateHelper(
        ql.QuoteHandle(FRAs[(n, m)]), n, m, settlementDays, calendar, ql.ModifiedFollowing, False, dayCounter
    )
    for n, m in FRAs.keys()
]
```

```python
dayCounter = ql.Actual360()
months = 3
futuresHelpers = [
    ql.FuturesRateHelper(
        ql.QuoteHandle(futures[d]),
        d,
        months,
        calendar,
        ql.ModifiedFollowing,
        True,
        dayCounter,
        ql.QuoteHandle(ql.SimpleQuote(0.0)),
    )
    for d in futures.keys()
]
```

The discount curve for the swaps will come from elsewhere. A real application would use some kind of risk-free curve; here we're using a flat one for convenience.

```python
discountTermStructure = ql.YieldTermStructureHandle(
    ql.FlatForward(settlementDate, 0.04, ql.Actual360()))
```

```python
settlementDays = 2
fixedLegFrequency = ql.Annual
fixedLegTenor = ql.Period(1, ql.Years)
fixedLegAdjustment = ql.Unadjusted
fixedLegDayCounter = ql.Thirty360()
floatingLegFrequency = ql.Quarterly
floatingLegTenor = ql.Period(3, ql.Months)
floatingLegAdjustment = ql.ModifiedFollowing
swapHelpers = [
    ql.SwapRateHelper(
        ql.QuoteHandle(swaps[(n, unit)]),
        ql.Period(n, unit),
        calendar,
        fixedLegFrequency,
        fixedLegAdjustment,
        fixedLegDayCounter,
        ql.Euribor3M(),
        ql.QuoteHandle(),
        ql.Period("0D"),
        discountTermStructure,
    )
    for n, unit in swaps.keys()
]
```

### Term structure construction

```python
forecastTermStructure = ql.RelinkableYieldTermStructureHandle()
```

```python
helpers = depositHelpers + futuresHelpers + swapHelpers[1:]
depoFuturesSwapCurve = ql.PiecewiseFlatForward(settlementDate, helpers, ql.Actual360())
```

```python
helpers = depositHelpers + fraHelpers + swapHelpers
depoFraSwapCurve = ql.PiecewiseFlatForward(settlementDate, helpers, ql.Actual360())
```

### Swap pricing

```python
swapEngine = ql.DiscountingSwapEngine(discountTermStructure)
```

```python
nominal = 1000000
length = 5
maturity = calendar.advance(settlementDate, length, ql.Years)
payFixed = True
```

```python
fixedLegFrequency = ql.Annual
fixedLegAdjustment = ql.Unadjusted
fixedLegDayCounter = ql.Thirty360()
fixedRate = 0.04
```

```python
floatingLegFrequency = ql.Quarterly
spread = 0.0
fixingDays = 2
index = ql.Euribor3M(forecastTermStructure)
floatingLegAdjustment = ql.ModifiedFollowing
floatingLegDayCounter = index.dayCounter()
```

```python
fixedSchedule = ql.Schedule(
    settlementDate,
    maturity,
    fixedLegTenor,
    calendar,
    fixedLegAdjustment,
    fixedLegAdjustment,
    ql.DateGeneration.Forward,
    False,
)
floatingSchedule = ql.Schedule(
    settlementDate,
    maturity,
    floatingLegTenor,
    calendar,
    floatingLegAdjustment,
    floatingLegAdjustment,
    ql.DateGeneration.Forward,
    False,
)
```

We'll build a 5-years swap starting spot...

```python
spot = ql.VanillaSwap(
    ql.VanillaSwap.Payer,
    nominal,
    fixedSchedule,
    fixedRate,
    fixedLegDayCounter,
    floatingSchedule,
    index,
    spread,
    floatingLegDayCounter,
)
spot.setPricingEngine(swapEngine)
```

...and one starting 1 year forward.

```python
forwardStart = calendar.advance(settlementDate, 1, ql.Years)
forwardEnd = calendar.advance(forwardStart, length, ql.Years)
fixedSchedule = ql.Schedule(
    forwardStart,
    forwardEnd,
    fixedLegTenor,
    calendar,
    fixedLegAdjustment,
    fixedLegAdjustment,
    ql.DateGeneration.Forward,
    False,
)
floatingSchedule = ql.Schedule(
    forwardStart,
    forwardEnd,
    floatingLegTenor,
    calendar,
    floatingLegAdjustment,
    floatingLegAdjustment,
    ql.DateGeneration.Forward,
    False,
)
```

```python
forward = ql.VanillaSwap(
    ql.VanillaSwap.Payer,
    nominal,
    fixedSchedule,
    fixedRate,
    fixedLegDayCounter,
    floatingSchedule,
    index,
    spread,
    floatingLegDayCounter,
)
forward.setPricingEngine(swapEngine)
```

We'll price them both on the bootstrapped curves.

This is the quoted 5-years market rate; we expect the fair rate of the spot swap to match it.


```python
print(swaps[(5, ql.Years)].value())
```

```python
def show(swap):
    print("NPV         = %.2f" % swap.NPV())
    print("Fair spread = %.4f %%" % (swap.fairSpread()*100))
    print("Fair rate   =  %.4f %%" % (swap.fairRate()*100))
```

These are the results for the 5-years spot swap on the deposit/futures/swap curve...

```python
forecastTermStructure.linkTo(depoFuturesSwapCurve)
show(spot)
```

...and these are on the deposit/fra/swap curve.

```python
forecastTermStructure.linkTo(depoFraSwapCurve)
show(spot)
```

The same goes for the 1-year forward swap, except for the fair rate not matching the spot rate.

```python
forecastTermStructure.linkTo(depoFuturesSwapCurve)
show(forward)
```

```python
forecastTermStructure.linkTo(depoFraSwapCurve)
show(forward)
```

Modifying the 5-years swap rate and repricing will change the results:

```python
swaps[(5, ql.Years)].setValue(0.046)
```

```python
forecastTermStructure.linkTo(depoFuturesSwapCurve)
```

```python
show(spot)
```

```python
show(forward)
```

```python
forecastTermStructure.linkTo(depoFraSwapCurve)
```

```python
show(spot)
```

```python
show(forward)
```
