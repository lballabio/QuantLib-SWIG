import sysconfig
import concurrent.futures
import time

import QuantLib as ql

if __name__ == "__main__":
    def pricing() -> float:
        todaysDate = ql.Date(15, ql.May, 2025)
        ql.Settings.instance().evaluationDate = todaysDate
        underlying = ql.SimpleQuote(7.0)
        dividendYield = ql.FlatForward(todaysDate, 0.05, ql.Actual365Fixed())
        volatility = ql.BlackConstantVol(todaysDate, ql.TARGET(), 0.10, ql.Actual365Fixed())
        riskFreeRate = ql.FlatForward(todaysDate, 0.05, ql.Actual365Fixed())

        process = ql.BlackScholesMertonProcess(
            ql.QuoteHandle(underlying),
            ql.YieldTermStructureHandle(dividendYield),
            ql.YieldTermStructureHandle(riskFreeRate),
            ql.BlackVolTermStructureHandle(volatility),
        )
        option = ql.VanillaOption(
            ql.PlainVanillaPayoff(ql.Option.Call, 8.0),
            ql.EuropeanExercise(ql.Date(17, ql.May, 2026))
        )
        engine = ql.FdBlackScholesVanillaEngine(process, 100, 8000)
        option.setPricingEngine(engine)
        return option.NPV()

    print(bool(sysconfig.get_config_var("Py_GIL_DISABLED")))

    for w in [1, 2, 4, 8, 16, 32, 64]:
        start = time.time()

        with concurrent.futures.ThreadPoolExecutor(max_workers=w) as executor:
            retVals = [executor.submit(pricing) for i in range(400)]

            for val in retVals:
                val.result()

        print(w, 400/(time.time()-start))