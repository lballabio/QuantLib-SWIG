"""
 Copyright (C) 2025 Hiroto Ogawa

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <https://www.quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
"""
import unittest
import math
import QuantLib as ql

LAG = 0
DC = ql.ActualActual(ql.ActualActual.ISDA)
CAL = ql.NullCalendar()

def flat_rate(rate):
    return ql.FlatForward(
        LAG, CAL, ql.makeQuoteHandle(rate), DC)

class PerpetualFuturesTest(unittest.TestCase):
    def test_perpetual_futures(self):
        """ Testing PerpetualFutures pricing by discounted cashflow method against analytic formulae. """
        val_date = ql.Date(20, 6, 2025)
        ql.Settings.instance().evaluationDate = val_date

        test_dicts = [
            # Discrete time
            {"payoff": ql.PerpetualFutures.Linear, "funding": ql.PerpetualFutures.FundingWithPreviousSpot,
             "freq": ql.Period(3, ql.Months), "spot": 10000., "domRate": 0.04, "btcYield": 0.02,
             "fundingRate": 0.01, "rateDiff": 0.005, "tol": 1.e-6},
            {"payoff": ql.PerpetualFutures.Linear, "funding": ql.PerpetualFutures.FundingWithCurrentSpot,
             "freq": ql.Period(3, ql.Months), "spot": 10000., "domRate": 0.04, "btcYield": 0.02,
             "fundingRate": 0.01, "rateDiff": 0.005, "tol": 1.e-6},
            {"payoff": ql.PerpetualFutures.Inverse, "funding": ql.PerpetualFutures.FundingWithPreviousSpot,
             "freq": ql.Period(3, ql.Months), "spot": 10000., "domRate": 0.04, "btcYield": 0.02,
             "fundingRate": 0.01, "rateDiff": 0.005, "tol": 1.e-6},
            {"payoff": ql.PerpetualFutures.Inverse, "funding": ql.PerpetualFutures.FundingWithCurrentSpot,
             "freq": ql.Period(3, ql.Months), "spot": 10000., "domRate": 0.04, "btcYield": 0.02,
             "fundingRate": 0.01, "rateDiff": 0.005, "tol": 1.e-6},
             # Continuous time
            {"payoff": ql.PerpetualFutures.Linear, "funding": ql.PerpetualFutures.FundingWithPreviousSpot,
             "freq": ql.Period(0, ql.Months), "spot": 10000., "domRate": 0.04, "btcYield": 0.02,
             "fundingRate": 0.2, "rateDiff": 0.005, "tol": 1.e-6},
            {"payoff": ql.PerpetualFutures.Inverse, "funding": ql.PerpetualFutures.FundingWithPreviousSpot,
             "freq": ql.Period(0, ql.Months), "spot": 10000., "domRate": 0.04, "btcYield": 0.02,
             "fundingRate": 0.2, "rateDiff": 0.005, "tol": 1.e-6},
        ]

        for d in test_dicts:
            pf = ql.PerpetualFutures(d["payoff"], d["funding"], d["freq"], CAL, DC)
            domYC = ql.YieldTermStructureHandle(flat_rate(d["domRate"]))
            btcYC = ql.YieldTermStructureHandle(flat_rate(d["btcYield"]))
            spot = ql.QuoteHandle(ql.SimpleQuote(d["spot"]))

            fundingTimes = [0.0]
            fundingRates = [d["fundingRate"]]
            ir_diffs = [d["rateDiff"]]
            engine = ql.DiscountingPerpetualFuturesEngine(domYC, btcYC, spot, fundingTimes, fundingRates, ir_diffs,
                                                      ql.DiscountingPerpetualFuturesEngine.PiecewiseConstant)
            # pricing
            pf.setPricingEngine(engine)
            npv = pf.NPV()

            # analytic price
            # for details, refer to
            # Perpetual Futures Pricing, Damien Ackerer, Julien Hugonnier, Urban Jermann, 2024
            # https://finance.wharton.upenn.edu/~jermann/AHJ-main-10.pdf
            period = d["freq"]
            length = float(period.length())
            unit = period.units()
            if unit == ql.Years:
                dt = length
            elif unit == ql.Months:
                dt = length / 12.
            elif unit == ql.Weeks:
                dt = length * 7. / 365.
            elif unit == ql.Days:
                dt = length / 365.
            elif unit == ql.Hours:
                dt = length / 365. / 24.
            elif unit == ql.Minutes:
                dt = length / 365. / 24. / 60.
            elif unit == ql.Seconds:
                dt = length / 365. / 24. / 60. / 60.
            else:
                raise RuntimeError("unknown funding frequency unit: " + str(unit))
            
            # Discrete time
            if length > 0:
                if d["payoff"] == ql.PerpetualFutures.Linear:
                    if d["funding"] == ql.PerpetualFutures.FundingWithPreviousSpot:
                        # Equation (12) in the above paper
                        expected = (
                            d["spot"] * (d["fundingRate"] - d["rateDiff"]) * math.exp(d["btcYield"] * dt) /
                            (math.exp(d["btcYield"] * dt) - math.exp(d["domRate"] * dt) + d["fundingRate"] * math.exp(d["btcYield"] * dt)))
                    elif d["funding"] == ql.PerpetualFutures.FundingWithCurrentSpot:
                        # at the end of "3 Perpetual futures pricing" in the above paper
                        expected = (
                            d["spot"] * (d["fundingRate"] - d["rateDiff"]) * math.exp(d["domRate"] * dt) /
                            (math.exp(d["btcYield"] * dt) - math.exp(d["domRate"] * dt) + d["fundingRate"] * math.exp(d["domRate"] * dt)))
                elif d["payoff"] == ql.PerpetualFutures.Inverse:
                    if d["funding"] == ql.PerpetualFutures.FundingWithPreviousSpot:
                        # "Proposition 2" in the above paper
                        expected = (
                            d["spot"] *
                            (math.exp(d["domRate"] * dt) - math.exp(d["btcYield"] * dt) + d["fundingRate"] * math.exp(d["domRate"] * dt)) /
                            (d["fundingRate"] - d["rateDiff"]) / math.exp(d["domRate"] * dt))
                    elif d["funding"] == ql.PerpetualFutures.FundingWithCurrentSpot:
                        expected = (
                        d["spot"] *
                        (math.exp(d["domRate"] * dt) - math.exp(d["btcYield"] * dt) + d["fundingRate"] * math.exp(d["btcYield"] * dt)) /
                        (d["fundingRate"] - d["rateDiff"]) / math.exp(d["btcYield"] * dt))
            else:
                # Continuous time
                if d["payoff"] == ql.PerpetualFutures.Linear:
                    # "Proposition 3" in the above paper
                    expected = d["spot"] * (d["fundingRate"] - d["rateDiff"]) / (d["btcYield"] - d["domRate"] + d["fundingRate"])
                elif d["payoff"] == ql.PerpetualFutures.Inverse:
                    # "Proposition 4" in the above paper
                    expected = d["spot"] * (d["domRate"] - d["btcYield"] + d["fundingRate"]) / (d["fundingRate"] - d["rateDiff"])
            failed_msg = f"Perpetual future price {npv} differs from analytic price {expected} with relative tolerance {d['tol']}\n"
            self.assertAlmostEqual(npv, expected, delta=d["tol"] * expected, msg=failed_msg)

if __name__ == '__main__':
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
