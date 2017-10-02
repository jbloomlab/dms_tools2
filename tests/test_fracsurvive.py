"""Tests ``dms2_fracsurvive``.

Written by Jesse Bloom."""


import sys
import os
import glob
import unittest
import subprocess
import numpy
import pandas


class TestFracSurvive(unittest.TestCase):
    """Runs ``dms2_fracsurvive`` on test data."""

    def setUp(self):
        """Gets files set up appropriately."""
        self.indir = os.path.abspath(os.path.join(
                os.path.dirname(__file__), './fracsurvive_input_files/'))
        self.expecteddir = os.path.join(self.indir, './expected_output/')
        self.testdir = os.path.abspath(os.path.join(
                os.path.dirname(__file__), './test_fracsurvive_files/'))
        if os.path.isdir(self.testdir):
            for f in glob.glob('{0}/*.csv'.format(self.testdir)):
                os.remove(f)
        else:
            os.mkdir(self.testdir)

        self.mock = os.path.join(self.indir, 'mock-1_codoncounts.csv')
        self.assertTrue(os.path.isfile(self.mock))
        self.sel = os.path.join(self.indir, 'H17L19-c1-1_codoncounts.csv')
        self.assertTrue(os.path.isfile(self.sel), "no {0}".format(self.sel))
        self.libfracsurvive = 0.021
        self.mincounts = [0, 10]


    def test_FracSurvive(self):
        """Runs ``dms_fracsurvive``."""
        for mincounts in self.mincounts:
            name = 'mincounts{0}'.format(mincounts)
            cmds = ['dms2_fracsurvive', 
                    '--mock', self.mock, 
                    '--sel', self.sel, 
                    '--name', name,
                    '--libfracsurvive', str(self.libfracsurvive),
                    '--mincount', str(mincounts),
                    '--outdir', self.testdir,
                    '--excludestop', 'yes',
                    ]
            sys.stderr.write('\nRunning the following command:\n{0}\n'
                    .format(' '.join(cmds)))
            subprocess.check_call(cmds)

            for (suffix, columntotest) in [
                    ('_mutfracsurvive.csv', 'mutfracsurvive'),
                    ('_sitefracsurvive.csv', 'avgfracsurvive'),
                    ('_sitefracsurvive.csv', 'maxfracsurvive'),
                    ]:
                f = os.path.join(self.testdir, name + suffix)
                self.assertTrue(os.path.isfile(f))
                df = pandas.read_csv(f)
                expected_f = os.path.join(self.expecteddir, os.path.basename(f))
                expected_df = pandas.read_csv(expected_f)
                expected_col = expected_df[columntotest]
                actual_col = df[columntotest]
                self.assertTrue(len(expected_col) == len(actual_col),
                        "length mismatch: {0} and {1}\n{2} and {3}".format(
                        f, expected_f, len(actual_col), len(expected_col)))
                self.assertTrue(numpy.allclose(
                        numpy.nan_to_num(actual_col), 
                        numpy.nan_to_num(expected_col)),
                        "mismatch: {0} and {1}\nmax diff: {2}"
                        .format(f, expected_f, 
                        (actual_col - expected_col).abs().max()))



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
