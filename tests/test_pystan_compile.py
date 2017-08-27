"""Tests compilation of ``pystan`` models."""


import unittest
import dms_tools2.prefs


class test_pyStanModels(unittest.TestCase):
    """Compiles ``pystan`` models without error control."""

    MODEL = dms_tools2.prefs.StanModelNoneErr

    def test_pystan_compilation(self):
        """Tests compilation of models."""
        m = self.MODEL(verbose=True)
        self.assertTrue(m is not None)


class test_pyStanModels_Same(test_pyStanModels):
    """Compiles ``pystan`` models 'same' error control."""

    MODEL = dms_tools2.prefs.StanModelSameErr

    def test_pystan_compilation(self):
        """Tests compilation of models."""
        m = self.MODEL(verbose=True)
        self.assertTrue(m is not None)


class test_pyStanModels_Different(test_pyStanModels):
    """Compiles ``pystan`` models 'different' error control."""

    MODEL = dms_tools2.prefs.StanModelDifferentErr

    def test_pystan_compilation(self):
        """Tests compilation of models."""
        m = self.MODEL(verbose=True)
        self.assertTrue(m is not None)

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
