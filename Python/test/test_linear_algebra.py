import QuantLib as ql
import unittest


class ArrayTest(unittest.TestCase):
    def test_math(self):
        a = ql.Array([1, 2, 3])
        b = ql.Array([4, 5, 6])
        self.assertEqual(-a, ql.Array([-1, -2, -3]))
        self.assertEqual(a + 1, ql.Array([2, 3, 4]))
        self.assertEqual(a + b, ql.Array([5, 7, 9]))
        self.assertEqual(b - 1, ql.Array([3, 4, 5]))
        self.assertEqual(b - a, ql.Array([3, 3, 3]))
        self.assertEqual(a * 2, ql.Array([2, 4, 6]))
        self.assertEqual(3 * a, ql.Array([3, 6, 9]))
        self.assertEqual(a * b, ql.Array([4, 10, 18]))
        self.assertEqual(
            a * ql.Matrix([[1, 2], [3, 4], [5, 6]]),
            ql.Array([22, 28])
        )
        self.assertEqual(b / 2, ql.Array([2, 2.5, 3]))
        self.assertEqual(b / a, ql.Array([4, 2.5, 2]))
        self.assertEqual(a @ b, 32)

    def test_compare(self):
        for v1 in ([1, 2], [1, 2, 3], [2, 3, 4]):
            for v2 in ([1, 2], [1, 2, 3], [2, 3, 4]):
                self.assertEqual(ql.Array(v1) == ql.Array(v2), v1 == v2)
                self.assertEqual(ql.Array(v1) != ql.Array(v2), v1 != v2)


if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)

