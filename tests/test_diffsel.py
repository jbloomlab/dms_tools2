"""Tests ``dms2_diffsel``.

Written by Jesse Bloom."""


import sys
import os
import glob
import unittest
import subprocess
import numpy
import pandas


class TestDiffSel(unittest.TestCase):
    """Runs ``dms2_diffsel`` on test data."""

    def setUp(self):
        """Gets files set up appropriately."""
        self.subdir = os.path.abspath(os.path.join(
                os.path.dirname(__file__), './diffsel_input_files/'))
        self.testdir = os.path.abspath(os.path.join(
                os.path.dirname(__file__), './test_diffsel_files/'))
        if os.path.isdir(self.testdir):
            for f in glob.glob('{0}/*.csv'.format(self.testdir)):
                os.remove(f)
        else:
            os.mkdir(self.testdir)
        self.mock = os.path.join(self.subdir, 'L1_mock_r1counts.csv')
        self.assertTrue(os.path.isfile(self.mock))
        self.concentrations = ['c1', 'c3']
        self.mincounts = [0, 10]
        self.selected = os.path.join(self.subdir,
                'L1_H17L19_{concentration}_r1counts.csv')
        self.expecteddir = os.path.join(self.subdir, './expected_output/')
        self.no_err_counts = os.path.join(self.subdir, 'no_err_counts.csv')
        self.err_counts = os.path.join(self.subdir, 'err_counts.csv')
        for c in self.concentrations:
            fname = self.selected.format(concentration=c)
            self.assertTrue(os.path.isfile(fname))
            for mincounts in self.mincounts:
                fprefix = os.path.join(self.expecteddir, 
                        '{0}-mincounts{1}'.format(c, mincounts))
                for fsuffix in ['_mutdiffsel.csv', '_sitediffsel.csv']:
                    fname = fprefix + fsuffix
                    self.assertTrue(os.path.isfile(fname), fname)


    def test_DiffSel(self):
        """Runs ``dms_diffsel``."""
        for c in self.concentrations:
            selected = self.selected.format(concentration=c)
            for mincounts in self.mincounts:
                name = '{0}-mincounts{1}'.format(c, mincounts)
                cmds = ['dms2_diffsel', 
                        '--mock', self.mock, 
                        '--sel', selected, 
                        '--name', name, 
                        '--mincount', str(mincounts),
                        '--pseudocount', '10',
                        '--outdir', self.testdir,
                        '--excludestop', 'yes',
                        ]
                sys.stderr.write('\nRunning the following command:\n{0}\n'
                        .format(' '.join(cmds)))
                subprocess.check_call(cmds)
                for (suffix, columntotest) in [
                        ('_mutdiffsel.csv', 'mutdiffsel'),
                        ('_sitediffsel.csv', 'abs_diffsel')
                        ]:
                    f = os.path.join(self.testdir, name + suffix)
                    self.assertTrue(os.path.isfile(f))
                    df = pandas.read_csv(f)
                    expected_f = os.path.join(self.expecteddir, name + suffix)
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
                            (actual_col - expected_col).max()))


    def test_ErrorCorrection_NoErrors(self):
        """Tests ``dms2_diffsel`` ``--err`` when no errors in control."""
        for c in self.concentrations:
            selected = self.selected.format(concentration=c)
            for mincounts in [0]:
                name = 'no-err-{0}-mincounts{1}'.format(c, mincounts)
                cmds = ['dms2_diffsel', 
                        '--mock', self.mock, 
                        '--sel', selected, 
                        '--err', self.no_err_counts,
                        '--name', name, 
                        '--mincount', str(mincounts),
                        '--pseudocount', '10',
                        '--outdir', self.testdir,
                        '--excludestop', 'yes',
                        ]
                sys.stderr.write('\nRunning the following command:\n{0}\n'
                        .format(' '.join(cmds)))
                subprocess.check_call(cmds)
                for (suffix, columntotest) in [
                        ('_mutdiffsel.csv', 'mutdiffsel'),
                        ('_sitediffsel.csv', 'abs_diffsel')
                        ]:
                    f = os.path.join(self.testdir, name + suffix)
                    self.assertTrue(os.path.isfile(f))
                    df = pandas.read_csv(f)
                    expected_f = os.path.join(self.expecteddir, 
                            name.replace('no-err-', '') + suffix)
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
                            (actual_col - expected_col).max()))


    def test_ErrorCorrection_WithErrors(self):
        """Tests ``dms2_diffsel`` ``--err`` wth non-zero error control."""
        for c in self.concentrations:
            selected = self.selected.format(concentration=c)
            for mincounts in [0]:
                name = 'errors-{0}-mincounts{1}'.format(c, mincounts)
                cmds = ['dms2_diffsel', 
                        '--mock', self.mock, 
                        '--sel', selected, 
                        '--err', self.err_counts,
                        '--name', name, 
                        '--mincount', str(mincounts),
                        '--pseudocount', '10',
                        '--outdir', self.testdir,
                        '--excludestop', 'yes',
                        ]
                sys.stderr.write('\nRunning the following command:\n{0}\n'
                        .format(' '.join(cmds)))
                subprocess.check_call(cmds)
                for (suffix, columntotest) in [
                        ('_mutdiffsel.csv', 'mutdiffsel'),
                        ('_sitediffsel.csv', 'abs_diffsel'),
                        ]:
                    f = os.path.join(self.testdir, name + suffix)
                    self.assertTrue(os.path.isfile(f))
                    df = pandas.read_csv(f)
                    expected_f = os.path.join(self.expecteddir, 
                            name + suffix)
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
                            (actual_col - expected_col).max()))



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
