# Example of option baskets
# Distributed under BSD License

from enthought.mayavi.scripts import mayavi2

mayavi2.standalone(globals())


import scipy
import QuantLib as ql
from enthought.tvtk.tools import mlab
from enthought.mayavi.sources.vtk_data_source import VTKDataSource
from enthought.mayavi.filters.warp_scalar import WarpScalar
from enthought.mayavi.modules.surface import Surface


spot = scipy.arange(10.0, 100.0, 5.0)
vol = scipy.arange(0.1, 1.0, 0.1)
riskfree = scipy.arange(0.0, 5.0, 1.0)

todaysDate = ql.Date(15, ql.May, 1998)
ql.Settings.instance().evaluationDate = todaysDate
settlementDate = ql.Date(17, ql.May, 1998)
riskFreeQuote = ql.SimpleQuote(0.05)
riskFreeRate = ql.FlatForward(settlementDate, ql.QuoteHandle(riskFreeQuote), ql.Actual365Fixed())

# option parameters
exercise1 = ql.AmericanExercise(settlementDate, ql.Date(17, ql.May, 1999))
exercise2 = ql.EuropeanExercise(settlementDate)
payoff = ql.PlainVanillaPayoff(ql.Option.Call, 40.0)

# market data
underlying = ql.SimpleQuote(36.0)
volatilityQuote = ql.SimpleQuote(0.05)
volatility = ql.BlackConstantVol(todaysDate, ql.QuoteHandle(volatilityQuote), ql.Actual365Fixed())
dividendYield = ql.FlatForward(settlementDate, 0.00, ql.Actual365Fixed())

# good to go
process = ql.BlackScholesMertonProcess(
    ql.QuoteHandle(underlying),
    ql.YieldTermStructureHandle(dividendYield),
    ql.YieldTermStructureHandle(riskFreeRate),
    ql.BlackVolTermStructureHandle(volatility),
)
option1 = ql.VanillaOption(process, payoff, exercise1)
option1.setPricingEngine(ql.BaroneAdesiWhaleyEngine())


def f(x, y):
    underlying.setValue(x)
    volatilityQuote.setValue(y)
    return option1.NPV()


def add_data(tvtk_data):
    """Add a TVTK data object `tvtk_data` to the mayavi pipleine.
    """
    d = VTKDataSource()
    d.data = tvtk_data
    mayavi.add_source(d)
    return d


def surf_regular(source):
    """Now visualize the data as done in mlab.
    """
    w = WarpScalar()
    source.add_child(w)
    s = Surface()
    w.add_child(s)


# 3D visualization of f:

if __name__ == "__main__":
    mayavi.new_scene()
    for r in riskfree:
        riskFreeQuote.setValue(r)
        s1 = mlab.SurfRegular(spot, vol, scipy.vectorize(f), scale=[1, 100, 1])
        s1.lut.alpha_range = (0.2, 0.2)
        d = add_data(s1.data)
        surf_regular(d)
