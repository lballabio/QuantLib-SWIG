import QuantLib as ql
import unittest


class SettingsTest(unittest.TestCase):
    def test_properties(self):
        settings = ql.Settings.instance()
        with ql.SavedSettings():
            for v in (ql.Date(1, 1, 2023), ql.Date(2, 2, 2023)):
                settings.evaluationDate = v
                self.assertEqual(settings.evaluationDate, v)
            for v in (True, False):
                settings.enforcesTodaysHistoricFixings = v
                self.assertEqual(settings.enforcesTodaysHistoricFixings, v)
            for v in (True, False):
                settings.includeReferenceDateEvents = v
                self.assertEqual(settings.includeReferenceDateEvents, v)
            for v in (True, False, None):
                settings.includeTodaysCashFlows = v
                self.assertEqual(settings.includeTodaysCashFlows, v)

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
