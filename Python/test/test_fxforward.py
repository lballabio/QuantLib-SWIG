import math
import pytest
import QuantLib as ql

# PyTest fixture for common setup (analogous to the C++ CommonVars struct)
@pytest.fixture
def common_vars():
    # Set evaluation date
    today = ql.Date(15, ql.March, 2024)
    ql.Settings.instance().evaluationDate = today

    # Define maturity date 6 months from today
    maturity_date = today + ql.Period(6, ql.Months)

    # Currencies
    usd = ql.USDCurrency()
    sgd = ql.SGDCurrency()

    # Flat yield term structures for USD and SGD discount curves
    # (using flat forward rates of 5.0% for USD and 3.5% for SGD)
    usd_curve = ql.FlatForward(today, 0.05, ql.Actual365Fixed())
    sgd_curve = ql.FlatForward(today, 0.035, ql.Actual365Fixed())
    usd_curve_handle = ql.RelinkableYieldTermStructureHandle(usd_curve)
    sgd_curve_handle = ql.RelinkableYieldTermStructureHandle(sgd_curve)

    # Spot FX quote: 1 USD = 1.35 SGD
    spot_quote = ql.SimpleQuote(1.35)
    spot_handle = ql.RelinkableQuoteHandle(spot_quote)

    # Tolerance for floating-point comparisons
    tolerance = 1.0e-6

    # Return all common variables in a dict for use by tests
    return {
        "today": today,
        "maturity_date": maturity_date,
        "usd": usd,
        "sgd": sgd,
        "usd_curve_handle": usd_curve_handle,
        "sgd_curve_handle": sgd_curve_handle,
        "spot_quote": spot_quote,
        "spot_handle": spot_handle,
        "tolerance": tolerance
    }

def test_fx_forward_construction(common_vars):
    """Testing FX forward construction with two nominal amounts."""
    vars = common_vars
    usd_nominal = 1_000_000.0
    sgd_nominal = 1_350_000.0

    # Create an FX forward (pay USD, receive SGD) with specified nominals
    fwd = ql.FxForward(usd_nominal, vars["usd"], sgd_nominal, vars["sgd"], 
                       vars["maturity_date"], True)  # paySourceCurrency=True means paying source (USD)

    # Check that attributes match the inputs
    assert fwd.sourceNominal() == pytest.approx(usd_nominal, rel=vars["tolerance"])
    assert fwd.targetNominal() == pytest.approx(sgd_nominal, rel=vars["tolerance"])
    assert fwd.sourceCurrency() == vars["usd"]
    assert fwd.targetCurrency() == vars["sgd"]
    assert fwd.maturityDate() == vars["maturity_date"]
    assert fwd.paySourceCurrency() is True  # paying the source currency (USD)
    assert fwd.isExpired() is False         # not expired since maturity is in the future

def test_fx_forward_construction_with_rate(common_vars):
    """Testing FX forward construction using a forward rate instead of explicit target nominal."""
    vars = common_vars
    nominal = 1_000_000.0
    forward_rate = 1.36  # forward FX rate (SGD per USD)

    # Create an FX forward (sell USD) by specifying a forward rate instead of target nominal.
    # paySourceCurrency=True -> paying USD, receiving SGD.
    fwd = ql.FxForward(nominal, vars["usd"], vars["sgd"], forward_rate, 
                       vars["maturity_date"], True)

    # The target nominal should be nominal * forward_rate (i.e., 1,000,000 * 1.36)
    expected_target_nominal = nominal * forward_rate
    assert fwd.sourceNominal() == pytest.approx(nominal, rel=vars["tolerance"])
    assert fwd.targetNominal() == pytest.approx(expected_target_nominal, rel=1e-6)
    assert fwd.sourceCurrency() == vars["usd"]
    assert fwd.targetCurrency() == vars["sgd"]
    # fwd.forwardRate() should equal the input forward_rate
    assert fwd.forwardRate() == pytest.approx(forward_rate, rel=1e-10)

def test_contracted_forward_rate(common_vars):
    """Testing that the FX forward reports the contracted forward rate correctly."""
    vars = common_vars
    usd_nominal = 1_000_000.0
    sgd_nominal = 1_350_000.0
    expected_rate = sgd_nominal / usd_nominal  # expected contracted rate = 1.35

    # FX forward constructed with explicit nominals
    fwd1 = ql.FxForward(usd_nominal, vars["usd"], sgd_nominal, vars["sgd"], 
                        vars["maturity_date"], True)
    assert fwd1.forwardRate() == pytest.approx(expected_rate, rel=1e-10)

    # FX forward constructed using a forward rate
    input_rate = 1.36
    fwd2 = ql.FxForward(usd_nominal, vars["usd"], vars["sgd"], input_rate, 
                        vars["maturity_date"], True)
    assert fwd2.forwardRate() == pytest.approx(input_rate, rel=1e-10)

    # Now attach a pricing engine and compare contracted vs fair forward rate
    engine = ql.DiscountingFxForwardEngine(vars["usd_curve_handle"], vars["sgd_curve_handle"], vars["spot_handle"])
    fwd1.setPricingEngine(engine)
    contracted_rate = fwd1.forwardRate()      # contractual forward rate (still 1.35 in fwd1)
    fair_rate = fwd1.fairForwardRate()        # fair forward rate given curves and spot
    # The contracted rate should generally differ from the fair forward rate (unless the deal was struck at fair value)
    assert contracted_rate != pytest.approx(fair_rate, rel=1e-12)

def test_fx_forward_expiry(common_vars):
    """Testing FX forward expiry flag for past maturity dates."""
    vars = common_vars
    past_date = vars["today"] - ql.Period(1, ql.Days)
    # Create a forward with maturity in the past (yesterday)
    expired_fwd = ql.FxForward(1_000_000.0, vars["usd"], 1_350_000.0, vars["sgd"], 
                               past_date, True)
    assert expired_fwd.isExpired() is True

def test_discounting_fx_forward_engine(common_vars):
    """Testing that the DiscountingFxForwardEngine produces an NPV and fair forward rate."""
    vars = common_vars
    usd_nominal = 1_000_000.0
    sgd_nominal = 1_350_000.0

    fwd = ql.FxForward(usd_nominal, vars["usd"], sgd_nominal, vars["sgd"], 
                       vars["maturity_date"], True)  # pay USD, receive SGD
    engine = ql.DiscountingFxForwardEngine(vars["usd_curve_handle"], vars["sgd_curve_handle"], vars["spot_handle"])
    fwd.setPricingEngine(engine)

    # Ensure NPV is calculated (should not be "null")
    npv = fwd.NPV()
    assert npv is not None and math.isfinite(npv)  # NPV should be a finite number
    # Ensure fair forward rate is positive
    fair_rate = fwd.fairForwardRate()
    assert fair_rate > 0.0

def test_fair_forward_rate(common_vars):
    """Testing fair forward rate calculation vs. manual expected value."""
    vars = common_vars
    usd_nominal = 1_000_000.0
    sgd_nominal = 1_350_000.0

    fwd = ql.FxForward(usd_nominal, vars["usd"], sgd_nominal, vars["sgd"], 
                       vars["maturity_date"], True)
    engine = ql.DiscountingFxForwardEngine(vars["usd_curve_handle"], vars["sgd_curve_handle"], vars["spot_handle"])
    fwd.setPricingEngine(engine)

    # Calculate the expected fair forward rate: Spot * (DF_sgd / DF_usd),
    # where DF_x = discount factor from settlement date to maturity for that currency.
    settlement_date = fwd.settlementDate()
    spot_fx = vars["spot_quote"].value()
    df_usd = vars["usd_curve_handle"].discount(vars["maturity_date"]) / vars["usd_curve_handle"].discount(settlement_date)
    df_sgd = vars["sgd_curve_handle"].discount(vars["maturity_date"]) / vars["sgd_curve_handle"].discount(settlement_date)
    expected_fair_rate = spot_fx * (df_sgd / df_usd)

    calculated_fair_rate = fwd.fairForwardRate()
    # The calculated fair rate should closely match the expected formula value
    assert calculated_fair_rate == pytest.approx(expected_fair_rate, rel=1e-4)

def test_at_the_money(common_vars):
    """Testing that an at-the-money FX forward (strike = fair forward rate) has zero NPV."""
    vars = common_vars
    spot_fx = vars["spot_quote"].value()
    usd_nominal = 1_000_000.0

    # Determine settlement date (for a default 2-day spot settlement forward).
    temp_fwd = ql.FxForward(usd_nominal, vars["usd"], usd_nominal, vars["sgd"], 
                             vars["maturity_date"], True)  # temporary to get settlementDate
    settlement_date = temp_fwd.settlementDate()

    # Compute discount factors from settlement to maturity for USD and SGD
    df_usd = vars["usd_curve_handle"].discount(vars["maturity_date"]) / vars["usd_curve_handle"].discount(settlement_date)
    df_sgd = vars["sgd_curve_handle"].discount(vars["maturity_date"]) / vars["sgd_curve_handle"].discount(settlement_date)
    fair_forward_rate = spot_fx * (df_sgd / df_usd)  # fair forward FX rate for reference

    # Compute the SGD nominal that makes the forward at-the-money: target_nominal = source_nominal * df_usd * spot / df_sgd
    sgd_nominal_atm = usd_nominal * df_usd * spot_fx / df_sgd

    # Create an at-the-money forward (pay USD, receive SGD at the fair forward rate)
    fwd_atm = ql.FxForward(usd_nominal, vars["usd"], sgd_nominal_atm, vars["sgd"], 
                           vars["maturity_date"], True)
    fwd_atm.setPricingEngine(ql.DiscountingFxForwardEngine(vars["usd_curve_handle"], vars["sgd_curve_handle"], vars["spot_handle"]))

    # NPV should be approximately zero for an ATM forward
    npv_atm = fwd_atm.NPV()
    assert abs(npv_atm) < 1.0e-4  # allow a tiny tolerance around zero

def test_position_direction(common_vars):
    """Testing that long vs. short positions yield opposite NPVs."""
    vars = common_vars
    usd_nominal = 1_000_000.0
    sgd_nominal = 1_350_000.0

    # Long USD position: pay SGD, receive USD (paySourceCurrency=False means paying target currency)
    fwd_long_usd = ql.FxForward(usd_nominal, vars["usd"], sgd_nominal, vars["sgd"], 
                                vars["maturity_date"], False)
    # Short USD position: pay USD, receive SGD (paySourceCurrency=True)
    fwd_short_usd = ql.FxForward(usd_nominal, vars["usd"], sgd_nominal, vars["sgd"], 
                                 vars["maturity_date"], True)

    engine = ql.DiscountingFxForwardEngine(vars["usd_curve_handle"], vars["sgd_curve_handle"], vars["spot_handle"])
    fwd_long_usd.setPricingEngine(engine)
    fwd_short_usd.setPricingEngine(engine)

    npv_long = fwd_long_usd.NPV()
    npv_short = fwd_short_usd.NPV()
    # The NPVs should be opposites (long vs short have equal and opposite value)
    assert npv_long == pytest.approx(-npv_short, rel=1e-4)

def test_ir_curve_sensitivity(common_vars):
    """Testing FX forward NPV sensitivity to interest rate curve shifts."""
    vars = common_vars
    usd_nominal = 1_000_000.0
    sgd_nominal = 1_350_000.0

    fwd = ql.FxForward(usd_nominal, vars["usd"], sgd_nominal, vars["sgd"], 
                       vars["maturity_date"], True)  # pay USD, receive SGD
    engine = ql.DiscountingFxForwardEngine(vars["usd_curve_handle"], vars["sgd_curve_handle"], vars["spot_handle"])
    fwd.setPricingEngine(engine)

    npv_base = fwd.NPV()
    # Bump USD interest rates up by 10 bps (0.05 -> 0.051)
    vars["usd_curve_handle"].linkTo(ql.FlatForward(vars["today"], 0.051, ql.Actual365Fixed()))
    npv_usd_up = fwd.NPV()

    # Reset USD curve and bump SGD rates up by 10 bps (0.035 -> 0.036)
    vars["usd_curve_handle"].linkTo(ql.FlatForward(vars["today"], 0.05, ql.Actual365Fixed()))
    vars["sgd_curve_handle"].linkTo(ql.FlatForward(vars["today"], 0.036, ql.Actual365Fixed()))
    npv_sgd_up = fwd.NPV()

    # When paying USD (source) and receiving SGD (target):
    # - Higher USD rates -> USD leg PV (payment) is discounted more -> less negative -> NPV increases.
    # - Higher SGD rates -> SGD leg PV (receipt) is discounted more -> less positive -> NPV decreases.
    assert npv_usd_up > npv_base
    assert npv_sgd_up < npv_base

def test_spot_fx_sensitivity(common_vars):
    """Testing FX forward NPV sensitivity to spot FX changes."""
    vars = common_vars
    usd_nominal = 1_000_000.0
    sgd_nominal = 1_350_000.0

    fwd = ql.FxForward(usd_nominal, vars["usd"], sgd_nominal, vars["sgd"], 
                       vars["maturity_date"], True)  # pay USD, receive SGD
    fwd.setPricingEngine(ql.DiscountingFxForwardEngine(vars["usd_curve_handle"], vars["sgd_curve_handle"], vars["spot_handle"]))

    npv_base = fwd.NPV()

    # Increase spot FX to 1.40 (USD strengthens against SGD)
    vars["spot_handle"].linkTo(ql.SimpleQuote(1.40))
    npv_spot_up = fwd.NPV()

    # Decrease spot FX to 1.30 (USD weakens against SGD)
    vars["spot_handle"].linkTo(ql.SimpleQuote(1.30))
    npv_spot_down = fwd.NPV()

    # Spot quote convention: spot = SGD per USD.
    # If USD strengthens (spot up), the SGD payoff is worth less in USD terms -> NPV should drop.
    # If USD weakens (spot down), the SGD payoff is worth more in USD terms -> NPV should rise.
    assert npv_spot_up < npv_base
    assert npv_spot_down > npv_base

@pytest.mark.skip(reason="QuantLib-Python does not expose `additionalResults` for instruments[2](https://quant.stackexchange.com/questions/27837/calculating-the-greeks-for-quantlib-python-swaptions).")
def test_additional_results(common_vars):
    """Testing FX forward additional results (internal engine outputs) - not directly available in Python."""
    vars = common_vars
    # In C++, fwd.additionalResults() provides keys like "spotFx", "sourceCurrencyDiscountFactor", etc.
    # Since this method isn't exposed in Python, this test is skipped. As a workaround, one could manually 
    # compute and verify those values if needed (see note below).
    pass

def test_settlement_days(common_vars):
    """Testing FX forward settlementDays attribute without a specific calendar."""
    vars = common_vars
    usd_nominal = 1_000_000.0
    sgd_nominal = 1_350_000.0

    # Overnight (O/N) forward: 0 settlement days
    fwd_on = ql.FxForward(usd_nominal, vars["usd"], sgd_nominal, vars["sgd"], 
                          vars["maturity_date"], True, 0)
    assert fwd_on.settlementDays() == 0
    assert fwd_on.settlementDate() == vars["today"]  # 0 days -> settlement today

    # Tom/Next (T/N) forward: 1 settlement day
    fwd_tn = ql.FxForward(usd_nominal, vars["usd"], sgd_nominal, vars["sgd"], 
                          vars["maturity_date"], True, 1)
    assert fwd_tn.settlementDays() == 1
    assert fwd_tn.settlementDate() == vars["today"] + ql.Period(1, ql.Days)

    # Spot forward (S/N): 2 settlement days (default)
    fwd_sn = ql.FxForward(usd_nominal, vars["usd"], sgd_nominal, vars["sgd"], 
                          vars["maturity_date"], True)  # default settlementDays=2
    assert fwd_sn.settlementDays() == 2
    assert fwd_sn.settlementDate() == vars["today"] + ql.Period(2, ql.Days)

def test_settlement_days_with_calendar(common_vars):
    """Testing FX forward settlement with a calendar (business day adjustment)."""
    vars = common_vars
    usd_nominal = 1_000_000.0
    sgd_nominal = 1_350_000.0

    cal = ql.TARGET()  # calendar that observes weekends
    # Choose an evaluation date that is a Friday (15 March 2024 is Friday)
    friday = ql.Date(15, ql.March, 2024)
    ql.Settings.instance().evaluationDate = friday

    # Forward with 2 settlement days using a calendar: settlement should skip the weekend
    fwd = ql.FxForward(usd_nominal, vars["usd"], sgd_nominal, vars["sgd"], 
                       vars["maturity_date"], True, 2, cal)
    expected_settlement = cal.advance(friday, ql.Period(2, ql.Days))
    assert fwd.settlementDate() == expected_settlement

    # Restore evaluation date
    ql.Settings.instance().evaluationDate = vars["today"]
