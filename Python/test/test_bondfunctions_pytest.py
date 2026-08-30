"""
 Copyright (C) 2009 Joseph Malicki
 Copyright (C) 2019 Prasad Somwanshi
 Copyright (C) 2023 Francois Botha
  Copyright (C) 2026 Chirag Desai

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

import QuantLib as ql
import pytest


@pytest.fixture
def bond_setup():
    """Fixture for bond test setup"""
    ql.Settings.instance().evaluationDate = ql.Date(2, 1, 2010)
    settlement_days = 3
    face_amount = 100.0
    redemption = 100.0
    issue_date = ql.Date(2, 1, 2008)
    maturity_date = ql.Date(2, 1, 2018)
    calendar = ql.UnitedStates(ql.UnitedStates.GovernmentBond)
    settlement_date = calendar.advance(
        ql.Settings.instance().evaluationDate, settlement_days, ql.Days)
    day_counter = ql.ActualActual(ql.ActualActual.Bond)
    sched = ql.Schedule(
        issue_date,
        maturity_date,
        ql.Period(ql.Semiannual),
        calendar,
        ql.Unadjusted,
        ql.Unadjusted,
        ql.DateGeneration.Backward,
        False,
    )
    coupons = [0.05]

    bond = ql.FixedRateBond(
        settlement_days,
        face_amount,
        sched,
        coupons,
        day_counter,
        ql.Following,
        redemption,
        issue_date,
    )

    flat_forward = ql.FlatForward(
        issue_date, coupons[0], day_counter, ql.Compounded, ql.Semiannual
    )
    term_structure_handle = ql.RelinkableYieldTermStructureHandle(flat_forward)
    bondEngine = ql.DiscountingBondEngine(term_structure_handle)
    bond.setPricingEngine(bondEngine)

    yield {
        'settlement_days': settlement_days,
        'face_amount': face_amount,
        'redemption': redemption,
        'issue_date': issue_date,
        'maturity_date': maturity_date,
        'calendar': calendar,
        'settlement_date': settlement_date,
        'day_counter': day_counter,
        'sched': sched,
        'coupons': coupons,
        'bond': bond,
        'flat_forward': flat_forward,
        'term_structure_handle': term_structure_handle
    }

    ql.Settings.instance().evaluationDate = ql.Date()
        'settlement_days': settlement_days,
        'face_amount': face_amount,
        'redemption': redemption,
        'issue_date': issue_date,
        'maturity_date': maturity_date,
        'calendar': calendar,
        'settlement_date': settlement_date,
        'day_counter': day_counter,
        'sched': sched,
        'coupons': coupons,
        'bond': bond,
        'flat_forward': flat_forward,
        'term_structure_handle': term_structure_handle
    }


def test_startDate(bond_setup):
    """ Testing BondFunctions startDate. """
    assert ql.BondFunctions.startDate(bond_setup['bond']) == bond_setup['issue_date']


def test_maturityDate(bond_setup):
    """ Testing BondFunctions maturityDate. """
    assert ql.BondFunctions.maturityDate(bond_setup['bond']) == bond_setup['maturity_date']


def test_isTradable(bond_setup):
    """ Testing BondFunctions isTradable. """
    assert ql.BondFunctions.isTradable(bond_setup['bond'], ql.Date(1, 6, 2010))
    assert not ql.BondFunctions.isTradable(bond_setup['bond'], ql.Date(1, 1, 2028))


def test_previousCashFlowDate(bond_setup):
    """ Testing BondFunctions previousCashFlowDate. """
    assert ql.BondFunctions.previousCashFlowDate(bond_setup['bond'], ql.Date(1, 6, 2010)) == ql.Date(4, 1, 2010)


def test_nextCashFlowDate(bond_setup):
    """ Testing BondFunctions nextCashFlowDate. """
    assert ql.BondFunctions.nextCashFlowDate(bond_setup['bond'], ql.Date(1, 6, 2010)) == ql.Date(2, 7, 2010)


def test_previousCashFlowAmount(bond_setup):
    """ Testing BondFunctions previousCashFlowAmount. """
    assert round(ql.BondFunctions.previousCashFlowAmount(bond_setup['bond'], ql.Date(1, 6, 2010)), 4) == 2.5


def test_nextCashFlowAmount(bond_setup):
    """ Testing BondFunctions nextCashFlowAmount. """
    assert round(ql.BondFunctions.nextCashFlowAmount(bond_setup['bond'], ql.Date(1, 6, 2010)), 4) == 2.5


def test_previousCouponRate(bond_setup):
    """ Testing BondFunctions previousCouponRate. """
    assert ql.BondFunctions.previousCouponRate(bond_setup['bond']) == 0.05


def test_nextCouponRate(bond_setup):
    """ Testing BondFunctions nextCouponRate. """
    assert ql.BondFunctions.nextCouponRate(bond_setup['bond']) == 0.05


def test_accrualStartDate(bond_setup):
    """ Testing BondFunctions accrualStartDate. """
    assert ql.BondFunctions.accrualStartDate(bond_setup['bond'], ql.Date(1, 6, 2010)) == ql.Date(2, 1, 2010)


def test_accrualEndDate(bond_setup):
    """ Testing BondFunctions accrualEndDate. """
    assert ql.BondFunctions.accrualEndDate(bond_setup['bond'], ql.Date(1, 6, 2010)) == ql.Date(2, 7, 2010)


def test_accrualPeriod(bond_setup):
    """ Testing BondFunctions accrualPeriod. """
    assert ql.BondFunctions.accrualPeriod(bond_setup['bond'], ql.Date(1, 6, 2010)) == 0.5


def test_accrualDays(bond_setup):
    """ Testing BondFunctions accrualDays. """
    assert ql.BondFunctions.accrualDays(bond_setup['bond'], ql.Date(1, 10, 2010)) == 184


def test_accruedPeriod(bond_setup):
    """ Testing BondFunctions accruedPeriod. """
    assert round(ql.BondFunctions.accruedPeriod(bond_setup['bond'], ql.Date(1, 6, 2010)), 8) == 0.41436464


def test_accruedDays(bond_setup):
    """ Testing BondFunctions accruedDays. """
    assert ql.BondFunctions.accruedDays(bond_setup['bond'], ql.Date(1, 6, 2010)) == 150


def test_accruedAmount(bond_setup):
    """ Testing BondFunctions accruedAmount. """
    assert round(ql.BondFunctions.accruedAmount(bond_setup['bond'], ql.Date(1, 6, 2010)), 8) == 2.0718232


def test_bps(bond_setup):
    """ Testing BondFunctions bps. """
    assert round(ql.BondFunctions.bps(bond_setup['bond'], bond_setup['flat_forward']), 8) == 0.06527501
    assert round(ql.BondFunctions.bps(bond_setup['bond'],
                                      ql.InterestRate(0.03, bond_setup['day_counter'], ql.Compounded, ql.Annual)), 8) == 0.07071951
    assert round(ql.BondFunctions.bps(bond_setup['bond'],
                                      0.03, bond_setup['day_counter'], ql.Compounded, ql.Annual), 8) == 0.07071951


def test_cleanPrice(bond_setup):
    """ Testing BondFunctions cleanPrice. """
    assert round(ql.BondFunctions.cleanPrice(bond_setup['bond'], bond_setup['flat_forward']), 4) == 99.9448
    assert round(ql.BondFunctions.cleanPrice(bond_setup['bond'],
                                             ql.InterestRate(0.03, bond_setup['day_counter'], ql.Compounded, ql.Annual)), 4) == 114.2806
    assert round(ql.BondFunctions.cleanPrice(bond_setup['bond'],
                                             0.03, bond_setup['day_counter'], ql.Compounded, ql.Annual), 4) == 114.2806


def test_atmRate(bond_setup):
    """ Testing BondFunctions atmRate. """
    assert round(ql.BondFunctions.atmRate(bond_setup['bond'], bond_setup['flat_forward'],
                                          bond_setup['settlement_date'],
                                          ql.BondPrice(99.94475138121548, ql.BondPrice.Clean)), 4) == 0.05


def test_bondYield(bond_setup):
    """ Testing BondFunctions bondYield. """
    assert round(ql.BondFunctions.bondYield(bond_setup['bond'],
                                            ql.BondPrice(110, ql.BondPrice.Clean),
                                            bond_setup['day_counter'], ql.Compounded, ql.Annual), 8) == 0.03582431


def test_duration(bond_setup):
    """ Testing BondFunctions duration. """
    assert round(ql.BondFunctions.duration(bond_setup['bond'],
                                           ql.InterestRate(0.03, bond_setup['day_counter'], ql.Compounded, ql.Annual)), 4) == 6.5835
    assert round(ql.BondFunctions.duration(bond_setup['bond'],
                                           0.03, bond_setup['day_counter'], ql.Compounded, ql.Annual), 4) == 6.5835


def test_convexity(bond_setup):
    """ Testing BondFunctions convexity. """
    assert round(ql.BondFunctions.convexity(bond_setup['bond'],
                                            ql.InterestRate(0.03, bond_setup['day_counter'], ql.Compounded, ql.Annual)), 4) == 54.3498
    assert round(ql.BondFunctions.convexity(bond_setup['bond'],
                                            0.03, bond_setup['day_counter'], ql.Compounded, ql.Annual), 4) == 54.3498


def test_basisPointValue(bond_setup):
    """ Testing BondFunctions basisPointValue. """
    assert round(ql.BondFunctions.basisPointValue(bond_setup['bond'],
                                                   0.03, bond_setup['day_counter'], ql.Compounded, ql.Annual,
                                                   bond_setup['settlement_date']), 8) == -0.07527271
    assert round(ql.BondFunctions.basisPointValue(bond_setup['bond'],
                                                   ql.InterestRate(0.03, bond_setup['day_counter'], ql.Compounded, ql.Annual),
                                                   bond_setup['settlement_date']), 8) == -0.07527271


def test_yieldValueBasisPoint(bond_setup):
    """ Testing BondFunctions yieldValueBasisPoint. """
    assert round(ql.BondFunctions.yieldValueBasisPoint(bond_setup['bond'],
                                                        ql.InterestRate(0.03, bond_setup['day_counter'], ql.Compounded, ql.Annual),
                                                        ql.Date(1, 9, 2010)), 10) == -1.44145e-05
    assert round(ql.BondFunctions.yieldValueBasisPoint(bond_setup['bond'],
                                                        0.03, bond_setup['day_counter'], ql.Compounded, ql.Annual,
                                                        ql.Date(1, 9, 2010)), 10) == -1.44145e-05


def test_zSpread(bond_setup):
    """ Testing BondFunctions zSpread. """
    assert round(ql.BondFunctions.zSpread(bond_setup['bond'],
                                          ql.BondPrice(87.5, ql.BondPrice.Clean),
                                          bond_setup['flat_forward'],
                                          bond_setup['day_counter'], ql.Compounded, ql.Annual), 8) == 0.02125053



if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    pytest.main([__file__, '-v'])
