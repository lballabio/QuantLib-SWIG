"""
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl

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
import unittest


class DateTest(unittest.TestCase):
    def setUp(self):
        pass

    def testArithmetics(self):
        "Testing date arithmetics"
        date = ql.Date_minDate()

        dold = date.dayOfMonth()
        mold = date.month()
        yold = date.year()

        while date < ql.Date_maxDate():
            date += 1

            d = date.dayOfMonth()
            m = date.month()
            y = date.year()

            # check if skipping any date
            if not (
                (d == dold + 1 and m == mold and y == yold)
                or (d == 1 and m == mold + 1 and y == yold)
                or (d == 1 and m == 1 and y == yold + 1)
            ):
                self.fail(
                    """
wrong day, month, year increment
    date: %(t)s
    day, month, year: %(d)d, %(m)d, %(y)d
    previous:         %(dold)d, %(mold)d, %(yold)d
                """
                    % locals()
                )
            dold = d
            mold = m
            yold = y

    def testHolidayList(self):
        """ Testing Calendar testHolidayList() method. """
        holidayLstFunction = ql.Calendar.holidayList(ql.Poland(), ql.Date(31, 12, 2014), ql.Date(3, 4, 2015), False)
        holidayLstManual = (ql.Date(1, 1, 2015), ql.Date(6, 1, 2015))
        # check if dates both from function and from manual imput are the same
        self.assertTrue(all([(a == b) for a, b in zip(holidayLstFunction, holidayLstManual)]))

    def tearDown(self):
        pass


if __name__ == "__main__":
    print("testing QuantLib " + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(DateTest, "test"))
    unittest.TextTestRunner(verbosity=2).run(suite)
