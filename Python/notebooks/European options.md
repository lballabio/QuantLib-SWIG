---
jupyter:
  jupytext:
    formats: ipynb,md,py:percent
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

# European options in QuantLib

Copyright (&copy;) 2004, 2005, 2006, 2007, 2020 StatPro Italia srl

This file is part of QuantLib, a free-software/open-source library for financial quantitative analysts and developers - https://www.quantlib.org/

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

### Option construction

```python
exercise = ql.EuropeanExercise(ql.Date(17, ql.May, 1999))
payoff = ql.PlainVanillaPayoff(ql.Option.Call, 8.0)
```

```python
option = ql.VanillaOption(payoff, exercise)
```

### Market data

```python
underlying = ql.SimpleQuote(7.0)
dividendYield = ql.FlatForward(todaysDate, 0.05, ql.Actual365Fixed())
volatility = ql.BlackConstantVol(todaysDate, ql.TARGET(), 0.10, ql.Actual365Fixed())
riskFreeRate = ql.FlatForward(todaysDate, 0.05, ql.Actual365Fixed())
```

### Processes and models

```python
process = ql.BlackScholesMertonProcess(
    ql.QuoteHandle(underlying),
    ql.YieldTermStructureHandle(dividendYield),
    ql.YieldTermStructureHandle(riskFreeRate),
    ql.BlackVolTermStructureHandle(volatility),
)
```

```python
hestonProcess = ql.HestonProcess(
    ql.YieldTermStructureHandle(riskFreeRate),
    ql.YieldTermStructureHandle(dividendYield),
    ql.QuoteHandle(underlying),
    0.1 * 0.1,
    1.0,
    0.1 * 0.1,
    0.0001,
    0.0,
)
hestonModel = ql.HestonModel(hestonProcess)
```

### Pricing


We'll collect tuples of method name, option value, estimated error, and discrepancy from the analytic formula.

```python
results = []
```

#### analytic formula

```python
option.setPricingEngine(ql.AnalyticEuropeanEngine(process))
value = option.NPV()
refValue = value

results.append(('Analytic', value, None, None))
```

#### Heston semi-analytic formula

```python
option.setPricingEngine(ql.AnalyticHestonEngine(hestonModel))
value = option.NPV()

results.append(('Heston analytic', value, None, abs(value - refValue)))
```

#### Heston COS method

```python
option.setPricingEngine(ql.COSHestonEngine(hestonModel))
value = option.NPV()

results.append(('Heston COS', value, None, abs(value - refValue)))
```

#### Integral method

```python
option.setPricingEngine(ql.IntegralEngine(process))
value = option.NPV()

results.append(('Integral', value, None, abs(value - refValue)))
```

#### Finite-difference method

```python
timeSteps = 801
gridPoints = 800
```

```python
option.setPricingEngine(ql.FdBlackScholesVanillaEngine(process, timeSteps, gridPoints))
value = option.NPV()

results.append(('Finite diff.', value, None, abs(value - refValue)))
```

#### Binomial method

```python
timeSteps = 801
```

```python
for tree in ["JR", "CRR", "EQP", "Trigeorgis", "Tian", "LR", "Joshi4"]:
    option.setPricingEngine(ql.BinomialVanillaEngine(process, tree, timeSteps))
    value = option.NPV()

    results.append((f'Binomial ({tree})', value, None, abs(value - refValue)))
```

#### Monte Carlo method

```python
option.setPricingEngine(ql.MCEuropeanEngine(process, "pseudorandom", timeSteps=1,
                                            requiredTolerance=0.02, seed=42))
value = option.NPV()

results.append(("Monte Carlo (pseudo-random)", value, option.errorEstimate(), abs(value - refValue)))
```

```python
option.setPricingEngine(ql.MCEuropeanEngine(process, "lowdiscrepancy", timeSteps=1,
                                            requiredSamples=32768))
value = option.NPV()

results.append(("Monte Carlo (low-discrepancy)", value, None, abs(value - refValue)))
```

#### Results

```python
df = pd.DataFrame(results,
                  columns=["Method", "Option value", "Error estimate", "Actual error"])
```

```python
df.style.hide_index()
```

The following displays the results when this is run as a Python script (in which case the cell above is not displayed).

```python
if 'get_ipython' not in globals():
    print(df)
```
