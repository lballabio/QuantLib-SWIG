import QuantLib as ql
import unittest


class DayCountersTest(unittest.TestCase):
    def test_bus252(self):
        """Test Business252 daycounter"""

        calendar = ql.UnitedStates(ql.UnitedStates.GovernmentBond)

        #
        # Check that SWIG signature for Business252 calendar allows to
        # pass custom calendar into the class constructor.  Old
        # QuantLib-SWIG versions allow only to create Business252
        # calendar with default constructor parameter (Brazil
        # calendar), and generate an exception when trying to pass a
        # custom calendar as a parameter. So we just check here that
        # no exception occurs.
        #
        ql.Business252(calendar)

    def test_hash(self):
        for dc1 in (ql.Actual360(), ql.Thirty365()):
            for dc2 in (ql.Actual360(), ql.Thirty365()):
                expected = dc1.name() == dc2.name()
                self.assertEqual(dc1 == dc2, expected)
                self.assertEqual(dc1 != dc2, not expected)
                self.assertEqual(hash(dc1) == hash(dc2), expected)


if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
