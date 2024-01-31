import QuantLib as ql
import operator
import unittest


class MoneyTest(unittest.TestCase):
    def test_order(self):
        ops = [operator.eq, operator.ne, operator.lt, operator.le, operator.gt, operator.ge]
        usd = lambda v: ql.Money(v, ql.USDCurrency())
        for m1 in (usd(1), usd(2)):
            for m2 in (usd(1), usd(2)):
                for op in ops:
                    self.assertEqual(op(m1,  m2), op(m1.value(), m2.value()))


if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
