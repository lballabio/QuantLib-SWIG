---
jupyter:
  jupytext:
    formats: md,py:percent
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

# American options

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

### Global parameters

```python
todaysDate = ql.Date(15, ql.May, 1998)
ql.Settings.instance().evaluationDate = todaysDate
```

```python
interactive = "get_ipython" in globals()
```

### Option construction

```python
exercise = ql.AmericanExercise(todaysDate, ql.Date(17, ql.May, 1999))
payoff = ql.PlainVanillaPayoff(ql.Option.Put, 40.0)
```

```python
option = ql.VanillaOption(payoff, exercise)
```

### Market data

```python
underlying = ql.SimpleQuote(36.0)
dividendYield = ql.FlatForward(todaysDate, 0.00, ql.Actual365Fixed())
volatility = ql.BlackConstantVol(todaysDate, ql.TARGET(), 0.20, ql.Actual365Fixed())
riskFreeRate = ql.FlatForward(todaysDate, 0.06, ql.Actual365Fixed())
```

```python
process = ql.BlackScholesMertonProcess(
    ql.QuoteHandle(underlying),
    ql.YieldTermStructureHandle(dividendYield),
    ql.YieldTermStructureHandle(riskFreeRate),
    ql.BlackVolTermStructureHandle(volatility),
)
```

### Pricing

We'll collect tuples of method name, option value, and estimated error from the analytic formula.

```python
results = []
```

#### Analytic approximations

```python
option.setPricingEngine(ql.BaroneAdesiWhaleyEngine(process))
results.append(("Barone-Adesi-Whaley", option.NPV()))
```

```python
option.setPricingEngine(ql.BjerksundStenslandEngine(process))
results.append(("Bjerksund-Stensland", option.NPV()))
```

#### Finite-difference method

```python
timeSteps = 801
gridPoints = 800
```

```python
option.setPricingEngine(ql.FdBlackScholesVanillaEngine(process, timeSteps, gridPoints))
results.append(("finite differences", option.NPV()))
```

#### Binomial method

```python
timeSteps = 801
```

```python
for tree in ["JR", "CRR", "EQP", "Trigeorgis", "Tian", "LR", "Joshi4"]:
    option.setPricingEngine(ql.BinomialVanillaEngine(process, tree, timeSteps))
    results.append(("Binomial (%s)" % tree, option.NPV()))
```

### Results

```python
df = pd.DataFrame(results, columns=["Method", "Option value"])
df.style.hide_index()
```

The following displays the results when this is run as a Python script (in which case the cell above is not displayed).

```python
if not interactive:
    print(df)
```
