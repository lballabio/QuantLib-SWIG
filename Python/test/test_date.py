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

import datetime
import operator
import QuantLib as ql
import unittest


class DateTest(unittest.TestCase):
    def testArithmetics(self):
        "Testing date arithmetics"
        today = ql.Date.todaysDate()
        date = today - ql.Period(30, ql.Years)
        end_date = today + ql.Period(30, ql.Years)

        dold = date.dayOfMonth()
        mold = date.month()
        yold = date.year()

        while date < end_date:
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

    def test_hash(self):
        for date1 in (ql.Date(1, 2, 2020), ql.Date(3, 4, 2022), ql.Date()):
            for date2 in (ql.Date(1, 2, 2020), ql.Date(3, 4, 2022), ql.Date()):
                expected = str(date1) == str(date2)
                self.assertEqual(date1 == date2, expected)
                self.assertEqual(date1 != date2, not expected)
                self.assertEqual(hash(date1) == hash(date2), expected)

    def test_order(self):
        ops = [operator.eq, operator.ne, operator.lt, operator.le, operator.gt, operator.ge]
        for date1 in (ql.Date(1, 2, 2020), ql.Date(3, 4, 2022), ql.Date()):
            for date2 in (ql.Date(1, 2, 2020), ql.Date(3, 4, 2022), ql.Date()):
                for op in ops:
                    self.assertEqual(op(date1,  date2),
                                     op(date1.serialNumber(), date2.serialNumber()))

    def testHolidayList(self):
        """ Testing Calendar testHolidayList() method. """
        holidayLstFunction = ql.Calendar.holidayList(ql.Poland(), ql.Date(31, 12, 2014), ql.Date(3, 4, 2015), False)
        holidayLstManual = (ql.Date(1, 1, 2015), ql.Date(6, 1, 2015))
        # check if dates both from function and from manual input are the same
        self.assertTrue(all([(a == b) for a, b in zip(holidayLstFunction, holidayLstManual)]))

    def testConversion(self):
        for m in range(1, 13):
            ql_date = ql.Date(m * 2, m, 2020)
            py_date = datetime.date(2020, m, m * 2)
            self.assertEqual(ql_date.to_date(), py_date)
            self.assertEqual(ql.Date.from_date(py_date), ql_date)
            # datetime works as well
            py_dt = datetime.datetime(2020, m, m * 2)
            self.assertEqual(ql.Date.from_date(py_dt), ql_date)

        with self.assertRaisesRegex(RuntimeError, "from_date requires a date"):
            ql.Date.from_date("2020-01-02")


class PeriodTest(unittest.TestCase):
    def test_hash(self):
        for per1 in (ql.Period("1D"), ql.Period("1W"), ql.Period("12M"), ql.Period()):
            for per2 in (ql.Period("1D"), ql.Period("1Y"), ql.Period()):
                expected = str(per1.normalized()) == str(per2.normalized())
                self.assertEqual(per1 == per2, expected)
                self.assertEqual(per1 != per2, not expected)
                self.assertEqual(hash(per1) == hash(per2), expected)

    def test_order(self):
        ops = [operator.eq, operator.ne, operator.lt, operator.le, operator.gt, operator.ge]
        for per1 in (ql.Period("1D"), ql.Period("1W"), ql.Period("12M")):
            for per2 in (ql.Period("1D"), ql.Period("1Y")):
                for op in ops:
                    self.assertEqual(op(per1,  per2), op(per2.frequency(), per1.frequency()))


if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
