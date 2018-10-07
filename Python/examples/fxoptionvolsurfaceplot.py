# -*- coding: utf-8 -*-
"""
 Copyright (C) 2018 Wojciech Ślusarski

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
"""

import QuantLib as ql
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.style.use('ggplot')

ref_date = ql.Date(3, 9,  2018)
dates = [ql.Date(3, 10, 2018),
         ql.Date(3, 12, 2018),
         ql.Date(3, 2, 2019)]
t_nodes = [ql.Actual365Fixed().dayCount(ref_date, t) / 365
           for t in dates]
atm_vols = [0.05 * t_nodes[0] ** 0.5,
            0.06 * t_nodes[1] ** 0.5,
            0.1  * t_nodes[2] ** 0.5]
rr25d = [0.02, 0.015, 0.01]
bf25d = [0.02, 0.01, 0.005]
day_counter = ql.Actual365Fixed()
cal = ql.JointCalendar(ql.TARGET(), ql.Poland())

fx = ql.QuoteHandle(ql.SimpleQuote(4.0))
fore = ql.RelinkableYieldTermStructureHandle(ql.FlatForward(ref_date,
                                                        0.005,
                                                        day_counter))
dom = ql.RelinkableYieldTermStructureHandle(ql.FlatForward(ref_date,
                                                        0.005,
                                                        day_counter))
vol_surface = ql.FxBlackVannaVolgaVolatilitySurface(ref_date, dates,
                                                    atm_vols, rr25d, bf25d,
                                                    day_counter, cal,
                                                    fx, dom, fore)
vol_surface.enableExtrapolation()

expiries = np.arange(20, 210) / 365
strikes = np.linspace(3.7, 4.4, 100)


interpolated_vols = np.array([vol_surface.blackVol(t, k)
                              for t in expiries
                                for k in strikes]).reshape(len(expiries),
                                                            len(strikes)).T
expiries, strikes = np.meshgrid(expiries, strikes)

fig = plt.figure()
ax = fig.gca(projection='3d')

surface = ax.plot_surface(expiries, strikes, interpolated_vols,
                          cmap='viridis')
ax.set_xlabel("Expiry")
ax.set_ylabel("Strike")
ax.set_zlabel("Volatility")

plt.show()


