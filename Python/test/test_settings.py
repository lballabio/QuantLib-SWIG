import QuantLib as ql
import unittest


class SettingsTest(unittest.TestCase):
    def test_saved_settings(self):
        settings = ql.Settings.instance()
        settings.evaluationDate = ql.Date(1, 1, 2023)
        with ql.SavedSettings():
            settings.evaluationDate = ql.Date(2, 2, 2023)
            self.assertEqual(settings.evaluationDate, ql.Date(2, 2, 2023))
        # Test that SavedSettings restores settings on exit.
        self.assertEqual(settings.evaluationDate, ql.Date(1, 1, 2023))


if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
