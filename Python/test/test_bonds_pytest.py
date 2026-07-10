"""
 Copyright (C) 2009 Joseph Malicki
 Copyright (C) 2019 Prasad Somwanshi
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

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def fixed_rate_bond():
    ql.Settings.instance().evaluationDate = ql.Date(2, 1, 2010)

    settlement_days = 3
    face_amount = 100.0
    redemption = 100.0
    issue_date = ql.Date(2, 1, 2008)
    maturity_date = ql.Date(2, 1, 2018)
    calendar = ql.UnitedStates(ql.UnitedStates.GovernmentBond)
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
    bond_engine = ql.DiscountingBondEngine(term_structure_handle)
    bond.setPricingEngine(bond_engine)

    yield {
        "bond": bond,
        "flat_forward": flat_forward,
        "settlement_days": settlement_days,
        "face_amount": face_amount,
        "redemption": redemption,
        "issue_date": issue_date,
        "maturity_date": maturity_date,
        "day_counter": day_counter,
        "coupons": coupons,
    }

    ql.Settings.instance().evaluationDate = ql.Date()


@pytest.fixture
def fixed_rate_bond_kwargs():
    settlement_days = 3
    face_amount = 100.0
    redemption = 100.0
    issue_date = ql.Date(2, 1, 2008)
    maturity_date = ql.Date(2, 1, 2018)
    calendar = ql.UnitedStates(ql.UnitedStates.GovernmentBond)
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

    return {
        "settlement_days": settlement_days,
        "face_amount": face_amount,
        "redemption": redemption,
        "issue_date": issue_date,
        "maturity_date": maturity_date,
        "day_counter": day_counter,
        "sched": sched,
        "coupons": coupons,
    }


# ---------------------------------------------------------------------------
# FixedRateBond tests
# ---------------------------------------------------------------------------

def test_frequency(fixed_rate_bond):
    """Testing FixedRateBond frequency() method."""
    assert fixed_rate_bond["bond"].frequency() == ql.Semiannual


def test_day_counter(fixed_rate_bond):
    """Testing FixedRateBond dayCounter() method."""
    assert fixed_rate_bond["bond"].dayCounter() == fixed_rate_bond["day_counter"]


def test_simple_inspectors(fixed_rate_bond):
    """Testing FixedRateBond simple inspectors."""
    b = fixed_rate_bond
    assert b["bond"].settlementDays() == b["settlement_days"]
    assert b["bond"].notional() == b["face_amount"]
    assert b["bond"].issueDate() == b["issue_date"]
    assert b["bond"].maturityDate() == b["maturity_date"]


# def test_settlement_value(fixed_rate_bond):
#     """Testing FixedRateBond settlement value."""
#     orig_date = ql.Settings.evaluationDate
#     ql.Settings.evaluationDate = fixed_rate_bond["issue_date"] + 1 * ql.Months
#     assert round(fixed_rate_bond["bond"].settlementValue(100.0), 4) == 102.3098
#     ql.Settings.evaluationDate = orig_date


def test_cash_flows(fixed_rate_bond):
    """Testing that the FixedRateBond gives the expected cash flows."""
    b = fixed_rate_bond
    assert [round(cf.amount(), 4) for cf in b["bond"].cashflows()] == (
        20 * [round(b["face_amount"] * b["coupons"][0] / 2, 4)] + [round(b["redemption"], 4)]
    )


def test_redemption(fixed_rate_bond):
    """Testing FixedRateBond redemption value and date."""
    b = fixed_rate_bond
    assert b["bond"].redemption().date() == b["maturity_date"]
    assert b["bond"].redemption().amount() == b["redemption"]


def test_redemptions(fixed_rate_bond):
    """Testing FixedRateBond redemptions."""
    b = fixed_rate_bond
    redemptions = b["bond"].redemptions()
    assert len(redemptions) == 1
    assert redemptions[0].date() == b["maturity_date"]
    assert redemptions[0].amount() == b["redemption"]


def test_notional(fixed_rate_bond):
    """Testing FixedRateBond notional values."""
    assert fixed_rate_bond["bond"].notional() == 100.0
    assert fixed_rate_bond["bond"].notionals() == (100.0, 0)


def test_next_coupon(fixed_rate_bond):
    """Testing FixedRateBond correct next coupon amount."""
    b = fixed_rate_bond
    assert b["bond"].nextCouponRate(b["issue_date"]) == 0.05


def test_prev_coupon(fixed_rate_bond):
    """Testing FixedRateBond correct previous coupon amount."""
    assert fixed_rate_bond["bond"].previousCouponRate() == 0.05


def test_clean_price(fixed_rate_bond):
    """Testing FixedRateBond clean price."""
    b = fixed_rate_bond
    bond, dc, issue = b["bond"], b["day_counter"], b["issue_date"]

    assert round(bond.cleanPrice(0.05, dc, ql.Compounded, ql.Semiannual, issue), 4) == 99.9964
    assert (
        round(bond.cleanPrice(0.05, dc, ql.Compounded, ql.Semiannual, issue + ql.Period(1, ql.Months)), 4)
        == 99.9921
    )
    assert (
        round(bond.cleanPrice(0.06, dc, ql.Compounded, ql.Semiannual, issue + ql.Period(1, ql.Months)), 4)
        == 92.5985
    )


def test_dirty_price(fixed_rate_bond):
    """Testing FixedRateBond dirty price."""
    b = fixed_rate_bond
    bond, dc, issue = b["bond"], b["day_counter"], b["issue_date"]

    assert round(bond.dirtyPrice(0.05, dc, ql.Compounded, ql.Semiannual, issue), 4) == 99.9964
    assert (
        round(bond.dirtyPrice(0.05, dc, ql.Compounded, ql.Semiannual, issue + ql.Period(1, ql.Months)), 4)
        == 100.4179
    )
    assert (
        round(bond.dirtyPrice(0.06, dc, ql.Compounded, ql.Semiannual, issue + ql.Period(1, ql.Months)), 4)
        == 93.0244
    )


def test_clean_price_from_z_spread(fixed_rate_bond):
    """Testing FixedRateBond clean price derived from Z-spread."""
    b = fixed_rate_bond
    assert (
        round(
            ql.cleanPriceFromZSpread(
                b["bond"],
                b["flat_forward"],
                0.01,
                b["day_counter"],
                ql.Compounded,
                ql.Semiannual,
                b["issue_date"] + ql.Period(1, ql.Months),
            ),
            4,
        )
        == 92.5926
    )


# ---------------------------------------------------------------------------
# FixedRateBond kwargs tests
# ---------------------------------------------------------------------------

def test_constructor_kwargs(fixed_rate_bond_kwargs):
    """Testing FixedRateBond constructor with keyword arguments."""
    k = fixed_rate_bond_kwargs
    bond = ql.FixedRateBond(
        settlementDays=k["settlement_days"],
        schedule=k["sched"],
        paymentDayCounter=k["day_counter"],
        issueDate=k["issue_date"],
        coupons=k["coupons"],
        faceAmount=k["face_amount"],
    )
    assert type(bond) is ql.FixedRateBond
    assert bond.dayCounter() == k["day_counter"]
    assert bond.settlementDays() == k["settlement_days"]
    assert bond.issueDate() == k["issue_date"]
    assert bond.maturityDate() == k["maturity_date"]
    assert bond.redemption().date() == k["maturity_date"]
    assert bond.redemption().amount() == k["redemption"]
    assert bond.notional(k["issue_date"]) == 100.0
    assert bond.notionals() == (100.0, 0)


# ---------------------------------------------------------------------------
# AmortizingFixedRateBond tests
# ---------------------------------------------------------------------------

def test_amortizing_fixed_rate_bond_interest_rates():
    # see AmortizingBondTest::testBrazilianAmortizingFixedRateBond
    # in the C++ test suite

    nominals = [
            1000       , 983.33300000, 966.66648898, 950.00019204,
            933.33338867, 916.66685434, 900.00001759, 883.33291726,
            866.66619177, 849.99933423, 833.33254728, 816.66589633,
            799.99937871, 783.33299165, 766.66601558, 749.99946306,
            733.33297499, 716.66651646, 699.99971995, 683.33272661,
            666.66624140, 649.99958536, 633.33294599, 616.66615618,
            599.99951997, 583.33273330, 566.66633377, 549.99954356,
            533.33290739, 516.66625403, 499.99963400, 483.33314619,
            466.66636930, 449.99984658, 433.33320226, 416.66634063,
            399.99968700, 383.33290004, 366.66635221, 349.99953317,
            333.33290539, 316.66626012, 299.99948151, 283.33271031,
            266.66594695, 249.99932526, 233.33262024, 216.66590450,
            199.99931312, 183.33277035, 166.66617153, 149.99955437,
            133.33295388, 116.66633464,  99.99973207,  83.33307672,
             66.66646137,  49.99984602,  33.33324734,  16.66662367
    ]

    expected_amortizations = [
            16.66700000, 16.66651102, 16.66629694, 16.66680337,
            16.66653432, 16.66683675, 16.66710033, 16.66672548,
            16.66685753, 16.66678695, 16.66665095, 16.66651761,
            16.66638706, 16.66697606, 16.66655251, 16.66648807,
            16.66645852, 16.66679651, 16.66699333, 16.66648520,
            16.66665604, 16.66663937, 16.66678981, 16.66663620,
            16.66678667, 16.66639952, 16.66679021, 16.66663617,
            16.66665336, 16.66662002, 16.66648780, 16.66677688,
            16.66652271, 16.66664432, 16.66686163, 16.66665363,
            16.66678696, 16.66654783, 16.66681904, 16.66662777,
            16.66664527, 16.66677860, 16.66677119, 16.66676335,
            16.66662168, 16.66670502, 16.66671573, 16.66659137,
            16.66654276, 16.66659882, 16.66661715, 16.66660049,
            16.66661924, 16.66660257, 16.66665534, 16.66661534,
            16.66661534, 16.66659867, 16.66662367, 16.66662367
    ]

    expected_coupons = [
            5.97950399, 4.85474255, 5.27619136, 5.18522454,
            5.33753111, 5.24221882, 4.91231709, 4.59116258,
            4.73037674, 4.63940686, 4.54843737, 3.81920094,
            4.78359948, 3.86733691, 4.38439657, 4.09359456,
            4.00262671, 4.28531030, 3.82068947, 3.55165259,
            3.46502778, 3.71720657, 3.62189368, 2.88388676,
            3.58769952, 2.72800044, 3.38838360, 3.00196900,
            2.91100034, 3.08940793, 2.59877059, 2.63809514,
            2.42551945, 2.45615766, 2.59111761, 1.94857222,
            2.28751141, 1.79268582, 2.19248291, 1.81913832,
            1.90625855, 1.89350716, 1.48110584, 1.62031828,
            1.38600825, 1.23425366, 1.39521333, 1.06968563,
            1.03950542, 1.00065409, 0.90968563, 0.81871706,
            0.79726493, 0.63678002, 0.57187676, 0.49829046,
            0.31177086, 0.27290565, 0.19062560, 0.08662552
    ]

    settlement_days = 0
    issue_date = ql.Date(2, ql.March, 2020)
    maturity_date = ql.Date(2, ql.March, 2025)

    schedule = ql.Schedule(
        issue_date,
        maturity_date,
        ql.Period(ql.Monthly),
        ql.Brazil(ql.Brazil.Settlement),
        ql.Unadjusted,
        ql.Unadjusted,
        ql.DateGeneration.Backward,
        False,
    )

    coupons = ql.FixedRateLeg(
        schedule,
        nominals=nominals,
        couponRates=[0.0675],
        dayCount=ql.Business252(ql.Brazil()),
        compounding=ql.Compounded,
        compoundingFrequency=ql.Annual,
        paymentAdjustment=ql.Following,
    )

    bond = ql.Bond(
        settlement_days,
        schedule.calendar(),
        issue_date,
        coupons,
    )

    cashflows = bond.cashflows()

    assert len(cashflows) == 2 * len(nominals)

    for i in range(len(nominals)):
        assert round(expected_coupons[i], 5) == round(cashflows[2 * i].amount(), 5)
        assert round(expected_amortizations[i], 5) == round(cashflows[2 * i + 1].amount(), 5)
