"""
 Copyright (C) 2026 Arihant Lodha

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
import QuantLib as ql


class ECBTest(unittest.TestCase):
    def test_is_ecb_date(self):
        "Testing ECB.isECBdate"
        # February 8, 2023 is a known ECB maintenance period start date
        ecb_date = ql.Date(8, ql.February, 2023)
        self.assertTrue(ql.ECB.isECBdate(ecb_date))
        # The day before is not
        self.assertFalse(ql.ECB.isECBdate(ecb_date - 1))

    def test_is_ecb_code(self):
        "Testing ECB.isECBcode"
        self.assertTrue(ql.ECB.isECBcode("FEB23"))
        self.assertFalse(ql.ECB.isECBcode("NOTECB"))

    def test_code_and_date_roundtrip(self):
        "Testing ECB code/date roundtrip"
        ecb_date = ql.Date(8, ql.February, 2023)
        code = ql.ECB.code(ecb_date)
        self.assertEqual(code, "FEB23")
        self.assertEqual(ql.ECB.date(code), ecb_date)

    def test_date_from_month_year(self):
        "Testing ECB.date(month, year)"
        # ECB.date(month, year) should return the same date as nextDate
        ecb_date = ql.Date(8, ql.February, 2023)
        self.assertEqual(ql.ECB.date(ql.February, 2023), ecb_date)

    def test_next_date(self):
        "Testing ECB.nextDate"
        before = ql.Date(1, ql.January, 2023)
        # First ECB date on or after January 1, 2023 is February 8, 2023
        self.assertEqual(ql.ECB.nextDate(before), ql.Date(8, ql.February, 2023))
        # Next after February 8, 2023 is March 22, 2023
        ecb_feb = ql.Date(8, ql.February, 2023)
        self.assertEqual(ql.ECB.nextDate(ecb_feb), ql.Date(22, ql.March, 2023))
        # nextDate also accepts an ECB code
        self.assertEqual(ql.ECB.nextDate("FEB23"), ql.Date(22, ql.March, 2023))

    def test_next_dates(self):
        "Testing ECB.nextDates"
        before = ql.Date(1, ql.January, 2023)
        dates = ql.ECB.nextDates(before)
        self.assertGreater(len(dates), 0)
        self.assertEqual(dates[0], ql.Date(8, ql.February, 2023))
        # All returned dates must pass isECBdate
        for d in dates:
            self.assertTrue(ql.ECB.isECBdate(d))

    def test_next_code(self):
        "Testing ECB.nextCode"
        before = ql.Date(1, ql.January, 2023)
        self.assertEqual(ql.ECB.nextCode(before), "FEB23")
        # nextCode also accepts an ECB code
        self.assertEqual(ql.ECB.nextCode("FEB23"), "MAR23")


if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
