"""Tests ``dms2_batch_diffsel``.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess
import random
import numpy
import pandas



class test_batch_diffsel(unittest.TestCase):
    """Runs ``dms2_batch_diffsel`` on test data."""

    def setUp(self):
        """Set up input data."""

        self.testdir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './test_batch_diffsel_files/')
        if not os.path.isdir(self.testdir):
            os.mkdir(self.testdir)

        self.indir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './diffsel_input_files/')

        self.summaryprefix = 'summary'

        self.names = ['c1', 'c3']

        self.expected = {}
        for name in self.names:
            for seltype in ['mutdiffsel', 'sitediffsel']:
                self.expected[(name, seltype)] = os.path.join(self.indir,
                        './expected_output/', 
                        'errors-{0}-mincounts0_{1}.csv'.format(name, seltype))
                assert os.path.isfile(self.expected[(name, seltype)])

        # define output files
        self.outfiles = [self.summaryprefix + suffix for suffix in
                ['.log'] +
                ['_' + avgtype + diffseltype + '.pdf'
                        for avgtype in ['mean', 'median']
                        for diffseltype in ['maxdiffsel', 'minmaxdiffsel',
                        'positivediffsel', 'totaldiffsel']]
                ] + \
                [self.summaryprefix + '_H17L19-' + suffix for suffix in
                [seltype + 'diffselcorr.pdf' for seltype in 
                        ['absolutesite', 'maxmut', 'mut', 'positivesite']] +
                [seltype + 'diffsel.csv' for seltype in
                        ['meanmut', 'meansite', 'medianmut', 'mediansite']]]

        self.outfiles = [os.path.join(self.testdir, f) for f in self.outfiles]

        for f in self.outfiles:
            if os.path.isfile(f):
                os.remove(f)


    def test_dms2_batch_diffsel(self):
        """Runs ``dms2_batch_diffsel`` on test data."""

        batchfile = os.path.join(self.testdir, 'batch.csv')
        df = pandas.DataFrame({
                'group':['H17L19', 'H17L19'],
                'name':self.names,
                'sel':['L1_H17L19_{0}_r1counts.csv'.format(n) 
                        for n in self.names],
                'mock':['L1_mock_r1counts.csv'] * 2,
                'err':['err_counts.csv'] * 2,
                })
        df.to_csv(batchfile, index=False)

        cmds = [
                'dms2_batch_diffsel',
                '--batchfile', batchfile,
                '--summaryprefix', self.summaryprefix,
                '--outdir', self.testdir,
                '--indir', self.indir,
                '--pseudocount', '10',
               ]
        sys.stderr.write('\nRunning the following command:\n{0}\n'.format(
                ' '.join(cmds)))
        subprocess.check_call(cmds)

        for f in self.outfiles:
            self.assertTrue(os.path.isfile(f), "Failed to create {0}".format(f))

        for name in self.names:
            for seltype in ['mutdiffsel', 'sitediffsel']:
                expected = pandas.read_csv(self.expected[(name, seltype)])
                actual = pandas.read_csv(os.path.join(self.testdir,
                        'H17L19-{0}_{1}.csv'.format(name, seltype)))
                self.assertTrue((expected.columns == actual.columns).all())
                for col in expected.columns:
                    if pandas.api.types.is_numeric_dtype(actual[col]):
                        self.assertTrue(numpy.allclose(expected[col],
                                actual[col], equal_nan=True))
                    else:
                        self.assertTrue((expected[col] == actual[col]).all())


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
