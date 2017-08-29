"""Tests ``dms2_batch_prefs``.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess
import random
import numpy
import pandas



class test_batch_prefs(unittest.TestCase):
    """Runs ``dms2_batch_bcsubamp`` on test data using ``--method bayesian``."""

    METHOD = 'bayesian'

    def setUp(self):
        """Set up input data."""

        self.testdir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './test_batch_prefs_files/')
        if not os.path.isdir(self.testdir):
            os.mkdir(self.testdir)

        self.indir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './prefs_input_files/')

        self.summaryprefix = 'summary-{0}'.format(self.METHOD)

        # define output files
        self.outfiles = [os.path.join(self.testdir, 
                self.summaryprefix + suffix) for suffix in
                ['.log', '_prefscorr.pdf', '_avgprefs.csv']]

        for f in self.outfiles:
            if os.path.isfile(f):
                os.remove(f)


    def test_dms2_batch_prefs(self):
        """Runs ``dms2_batch_prefs`` on test data."""

        batchfile = os.path.join(self.testdir, 'batch.csv')
        replicates = [1, 2]
        df = pandas.DataFrame({
                'name':['{0}-replicate-{1}'.format(self.METHOD, r)
                        for r in replicates],
                'pre':['mutDNA-{0}_codoncounts.csv'.format(r)
                        for r in replicates],
                'post':['mutvirus-{0}_codoncounts.csv'.format(r)
                        for r in replicates],
                })
        df.to_csv(batchfile, index=False)

        cmds = [
                'dms2_batch_prefs',
                '--batchfile', batchfile,
                '--summaryprefix', self.summaryprefix,
                '--outdir', self.testdir,
                '--indir', self.indir,
                '--method', self.METHOD,
                '--ncpus', '1',
               ]
        sys.stderr.write('\nRunning the following command:\n{0}\n'.format(
                ' '.join(cmds)))
        subprocess.check_call(cmds)

        for f in self.outfiles:
            self.assertTrue(os.path.isfile(f), "Failed to create {0}".format(f))

        for fbase in ['{0}-replicate-1_prefs.csv'.format(self.METHOD),
                      '{0}-replicate-2_prefs.csv'.format(self.METHOD),
                      '{0}_avgprefs.csv'.format(self.summaryprefix)]:
            actual = os.path.join(self.testdir, fbase)
            self.assertTrue(os.path.isfile(actual), "No {0}".format(actual))
            expected = os.path.join(self.indir, fbase)
            actual_df = pandas.read_csv(actual, index_col='site')
            expected_df = pandas.read_csv(expected, index_col='site')
            for c in actual_df.columns:
                self.assertTrue(numpy.allclose(actual_df[c].values,
                        expected_df[c].values, atol=0.01),
                        '{0}, {1}, {2}'.format(self.METHOD, c,
                        numpy.abs(actual_df[c].values - 
                        expected_df[c].values).max()))


class test_batch_prefs_ratio(test_batch_prefs):
    """Runs ``dms2_batch_bcsubamp`` on test data using ``--method ratio``."""

    METHOD = 'ratio'


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
