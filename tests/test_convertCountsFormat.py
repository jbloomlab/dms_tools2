"""Tests `dms_tools2.utils.convertCountsFormat` function."""


import os
import unittest
import numpy
import pandas
from dms_tools2 import CODONS
import dms_tools2.utils


class test_convertCountsFormat(unittest.TestCase):
    """Tests `dms_tools2.utilsConvertCountsFormat`."""

    def test_convertCountsFormat(self):
        """Tests `dms_tools2.utilsConvertCountsFormat`."""
        testdir = os.path.join(os.path.dirname(__file__),
                'convertCountsFormat_input_files')
        outdir = os.path.join(os.path.dirname(__file__),
                'test_convertCountsFormat_files')
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        oldfile = os.path.join(testdir, 'dms_tools_format.txt')
        expectedfile = os.path.join(testdir, 'dms_tools2_format.csv')
        newfile = os.path.join(outdir, 'dms_tools2_format.csv')
        if os.path.isfile(newfile):
            os.remove(newfile)

        self.assertFalse(os.path.isfile(newfile))
        dms_tools2.utils.convertCountsFormat(oldfile, newfile, CODONS)
        self.assertTrue(os.path.isfile(newfile))

        new = pandas.read_csv(newfile)
        expected = pandas.read_csv(expectedfile)
        self.assertTrue(all(new.columns == expected.columns))
        for col in new.columns:
            self.assertTrue(all(new[col] == expected[col]))



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
