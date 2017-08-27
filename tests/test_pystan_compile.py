"""Tests compilation of ``pystan`` models."""


import unittest
import dms_tools2.prefs


class test_pyStanModels(unittest.TestCase):
    """Tests compilation of ``pystan`` models."""

    def test_pystan_compilation(self):
        """Tests compilation of models."""

        m = dms_tools2.prefs.StanModelNoneErr(verbose=True)
    
        m = dms_tools2.prefs.StanModelSameErr(verbose=True)

        m = dms_tools2.prefs.StanModelDifferentErr(verbose=True)


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
