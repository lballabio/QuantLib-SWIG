# Copyright (C) 2026 QuantLib contributors
# Copyright (C) 2026 Kyrylo Protsenko

# This file is part of QuantLib, a free-software/open-source library
# for financial quantitative analysts and developers - http://quantlib.org/

import unittest

import QuantLib as ql


class MtMCrossCurrencyBasisSwapTest(unittest.TestCase):

    def setUp(self):
        self.today = ql.Date(15, ql.January, 2025)
        ql.Settings.instance().evaluationDate = self.today
        self.calendar = ql.TARGET()
        self.curve = ql.YieldTermStructureHandle(
            ql.FlatForward(self.today, 0.02, ql.Actual365Fixed()))
        self.eur_index = ql.Euribor3M(self.curve)
        self.usd_index = ql.USDLibor(ql.Period(3, ql.Months), self.curve)

    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()

    def testFxResetCashFlowsAndPricer(self):
        """Testing FX-reset values, cash flows and their discounting pricer..."""
        convention = ql.FxResetConvention(2, self.calendar)
        value_date = ql.Date(17, ql.January, 2025)
        reset = convention.reset(value_date)

        self.assertEqual(reset.fixingDate(), self.today)
        self.assertEqual(reset.valueDate(), value_date)
        self.assertEqual(convention.valueDate(reset.fixingDate()), value_date)
        self.assertEqual(convention.reset(value_date).fixingDate(), self.today)
        self.assertEqual(convention.reset(value_date).valueDate(), value_date)

        pricer = ql.DiscountingFxResetPricer(
            ql.EURCurrency(), ql.USDCurrency(), self.curve, self.curve,
            ql.makeQuoteHandle(1.10), True)
        self.assertAlmostEqual(pricer.fxRate(reset), 1.10)

        exchange = ql.FxResetNotionalExchange(
            value_date, 100.0, None, reset)
        self.assertFalse(exchange.hasPreviousReset())
        self.assertTrue(exchange.hasCurrentReset())
        self.assertEqual(exchange.currentReset().fixingDate(), self.today)
        self.assertEqual(exchange.currentReset().valueDate(), value_date)

        exchange.setFxResetPricer(pricer)
        self.assertAlmostEqual(exchange.amount(), -110.0)

    def testSwapAndRateHelperExposeFxResetConventions(self):
        """Testing MtM swap and rate-helper FX-reset conventions..."""
        start = ql.Date(17, ql.January, 2025)
        end = ql.Date(17, ql.January, 2026)
        schedule = ql.Schedule(
            start, end, ql.Period(3, ql.Months), self.calendar,
            ql.Following, ql.Following, ql.DateGeneration.Forward, False)
        convention = ql.FxResetConvention(2, self.calendar)

        swap = ql.MtMCrossCurrencyBasisSwap(
            ql.MtMCrossCurrencyBasisSwap.Type_PayFxBaseCurrency,
            100.0, ql.EURCurrency(), schedule, self.eur_index, 0.0, 1.0,
            110.0, ql.USDCurrency(), schedule, self.usd_index, 0.0, 1.0,
            True, convention, 2, 3, ql.ModifiedFollowing, ql.Preceding,
            useIndexedCoupons=False)

        self.assertEqual(swap.fxResetConvention().fixingDays(), 2)
        self.assertEqual(swap.fxBasePaymentConvention(), ql.ModifiedFollowing)
        self.assertEqual(swap.fxQuotePaymentConvention(), ql.Preceding)
        self.assertEqual(swap.resettingLegIndex(), 0)
        self.assertEqual(
            len(swap.resettingLeg()), 2 * (len(schedule) - 1) + 1)

        coupons = [ql.as_fx_reset_coupon(cf) for cf in swap.resettingLeg()]
        coupons = [coupon for coupon in coupons if coupon]
        self.assertTrue(coupons)
        exchanges = [
            ql.as_fx_reset_notional_exchange(cf)
            for cf in swap.resettingLeg()]
        self.assertTrue(any(exchange for exchange in exchanges))
        first_coupon = coupons[0]
        self.assertEqual(
            first_coupon.fxResetValueDate(), first_coupon.accrualStartDate())
        self.assertEqual(
            first_coupon.fxResetDate(),
            self.calendar.advance(
                first_coupon.fxResetValueDate(), -2, ql.Days))

        swap.setPricingEngine(ql.DiscountingMtMCrossCurrencyBasisSwapEngine(
            ql.USDCurrency(), self.curve, ql.EURCurrency(), self.curve,
            ql.makeQuoteHandle(1.10)))
        reset_rates = swap.fxResetRates()
        reset_notionals = swap.fxResetNotionals()
        self.assertEqual(len(reset_rates), len(coupons))
        self.assertEqual(len(reset_notionals), len(coupons))
        self.assertEqual(swap.legCurrency(0), ql.EURCurrency())
        self.assertEqual(swap.legCurrency(1), ql.USDCurrency())
        self.assertAlmostEqual(swap.legNPV(0), 1.10 * swap.inCcyLegNPV(0))
        self.assertAlmostEqual(swap.legNPV(1), swap.inCcyLegNPV(1))
        for rate, notional in zip(reset_rates, reset_notionals):
            self.assertAlmostEqual(rate, 1.0 / 1.10)
            self.assertAlmostEqual(notional, 100.0)

        pricer = ql.DiscountingFxResetPricer(
            ql.USDCurrency(), ql.EURCurrency(), self.curve, self.curve,
            ql.makeQuoteHandle(1.0), True)
        ql.setFxResetPricer(swap.resettingLeg(), pricer)
        self.assertTrue(first_coupon.fxResetPricer())

        helper = ql.MtMCrossCurrencyBasisSwapRateHelper(
            ql.makeQuoteHandle(-0.001), ql.Period(1, ql.Years), 2,
            self.calendar, ql.Following, False, self.eur_index,
            self.usd_index, self.curve, False, True, True,
            ql.NoFrequency, 0, ql.NoFrequency, 2, self.calendar,
            False)
        self.assertEqual(helper.fxResetConvention().fixingDays(), 2)
        self.assertEqual(helper.swap().fxResetConvention().fixingDays(), 2)
        self.assertEqual(helper.swap().legCurrency(0), ql.EURCurrency())
        self.assertEqual(helper.swap().legCurrency(1), ql.USDCurrency())
        self.assertTrue(hasattr(helper.swap(), 'inCcyLegNPV'))

    def testIndexedCouponArguments(self):
        """Testing indexed-coupon arguments for cross-currency swaps..."""
        start = ql.Date(17, ql.January, 2025)
        end = ql.Date(17, ql.January, 2026)
        schedule = ql.Schedule(
            start, end, ql.Period(3, ql.Months), self.calendar,
            ql.Following, ql.Following, ql.DateGeneration.Forward, False)

        fixed_vs_floating_swap = (
            ql.ConstNotionalCrossCurrencyFixedVsFloatingSwap(
                ql.Swap.Payer,
                100.0, ql.EURCurrency(), schedule, 0.02,
                ql.Actual365Fixed(), ql.Following, 0, self.calendar,
                110.0, ql.USDCurrency(), schedule, self.usd_index, 0.0,
                ql.Following, 0, self.calendar,
                useIndexedCoupons=True))
        basis_swap = ql.ConstNotionalCrossCurrencyBasisSwap(
            100.0, ql.EURCurrency(), schedule, self.eur_index, 0.0, 1.0,
            110.0, ql.USDCurrency(), schedule, self.usd_index, 0.0, 1.0,
            useIndexedCoupons=False)
        mtm_swap = ql.MtMCrossCurrencyBasisSwap(
            ql.MtMCrossCurrencyBasisSwap.Type_PayFxBaseCurrency,
            100.0, ql.EURCurrency(), schedule, self.eur_index, 0.0, 1.0,
            110.0, ql.USDCurrency(), schedule, self.usd_index, 0.0, 1.0,
            False, useIndexedCoupons=True)

        quote = ql.makeQuoteHandle(-0.001)
        tenor = ql.Period(1, ql.Years)
        fixed_vs_floating_helper = (
            ql.ConstNotionalCrossCurrencySwapRateHelper(
                ql.makeQuoteHandle(0.02), tenor, 2, self.calendar,
                ql.Following, False, ql.Annual, ql.Actual365Fixed(),
                self.usd_index, self.curve, False, 0, True))
        basis_helper = ql.ConstNotionalCrossCurrencyBasisSwapRateHelper(
            quote, tenor, 2, self.calendar, ql.Following, False,
            self.eur_index, self.usd_index, self.curve, False, True,
            ql.Semiannual, 0, ql.Semiannual, False)
        mtm_helper = ql.MtMCrossCurrencyBasisSwapRateHelper(
            quote, tenor, 2, self.calendar, ql.Following, False,
            self.eur_index, self.usd_index, self.curve, False, True, False,
            ql.Semiannual, 0, ql.Semiannual, 0, self.calendar, True)

        self.assertTrue(fixed_vs_floating_swap)
        self.assertTrue(basis_swap)
        self.assertTrue(mtm_swap)
        self.assertTrue(fixed_vs_floating_helper.swap())
        self.assertTrue(basis_helper.swap())
        self.assertTrue(mtm_helper.swap())


if __name__ == '__main__':
    unittest.main()
