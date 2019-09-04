"""Tests ``dms2_batch_fracsurvive``.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess
import random
import numpy
import pandas


class test_batch_fracsurvive(unittest.TestCase):
    """Runs ``dms2_batch_fracsurvive`` on test data."""

    def setUp(self):
        """Set up input data."""

        self.testdir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './test_batch_fracsurvive_files/')
        if not os.path.isdir(self.testdir):
            os.mkdir(self.testdir)

        self.indir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './fracsurvive_input_files/')
        self.expectdir = os.path.join(self.indir, 'expected_output')

        self.summaryprefix = 'summary'

        # define output files
        self.outfiles = [self.summaryprefix + suffix for suffix in
                ['.log'] +
                ['_' + avgtype + fracsurvivetype + '.pdf'
                        for avgtype in ['mean', 'median']
                        for fracsurvivetype in ['avgfracsurvive', 
                        'maxfracsurvive']]
                ] + \
                [self.summaryprefix + '_' + antibody + '-' + seltype +
                        'fracsurvivecorr.pdf'
                    for seltype in ['avg', 'max', 'mut'] 
                    for antibody in ['H17L19-c1', 'H17L19-c3']] + \
                [self.summaryprefix + '_' + antibody + '-' + seltype +
                        'fracsurvive.csv' for seltype in
                        ['meanmut', 'meansite', 'medianmut', 'mediansite']
                        for antibody in ['H17L19-c1', 'H17L19-c3']]

        self.outfiles = [os.path.join(self.testdir, f) for f in self.outfiles]

        for f in self.outfiles:
            if os.path.isfile(f):
                os.remove(f)


    def test_dms2_batch_fracsurvive(self):
        """Runs ``dms2_batch_fracsurvive`` on test data."""

        batchfile = os.path.join(self.testdir, 'batch.csv')
        df = pandas.DataFrame.from_records(
                [(antibody,
                  'replicate-{0}'.format(r),
                  '{0}-{1}'.format(antibody, r),
                  'mock-{0}'.format(r),
                  'err-ctrl',
                  str(gamma),
                  antibodylabel
                  ) 
                 for r in [1, 2, 3]
                 for (antibody, gamma, antibodylabel) in reversed([
                    ('H17L19-c1', 0.021, r'H17L19 1 $\mu$g/ml'), 
                    ('H17L19-c3', 0.003, r'H17L19 10 $\mu$g/ml')])
                ],
                columns=['group', 'name', 'sel', 'mock', 'err', 
                    'libfracsurvive', 'grouplabel']
                )
        df.to_csv(batchfile, index=False)

        cmds = [
                'dms2_batch_fracsurvive',
                '--batchfile', batchfile,
                '--summaryprefix', self.summaryprefix,
                '--outdir', self.testdir,
                '--indir', self.indir,
               ]
        sys.stderr.write('\nRunning the following command:\n{0}\n'.format(
                ' '.join(cmds)))
        subprocess.check_call(cmds)

        for f in self.outfiles:
            self.assertTrue(os.path.isfile(f), "Failed to create {0}".format(f))

        for seltype in ['mutfracsurvive', 'sitefracsurvive']:
            for antibody in ['H17L19-c1', 'H17L19-c3']:
                for avgtype in ['mean', 'median']:
                    base_f = '{0}_{1}-{2}{3}.csv'.format(
                            self.summaryprefix, antibody,
                            avgtype, seltype)
                    expected = pandas.read_csv(os.path.join(
                            self.expectdir, base_f))
                    actual = pandas.read_csv(os.path.join(
                            self.testdir, base_f))
                    self.assertTrue((expected.columns == actual.columns).all())
                    for col in actual.columns:
                        if pandas.api.types.is_numeric_dtype(actual[col]):
                            self.assertTrue(numpy.allclose(expected[col],
                                    actual[col], equal_nan=True))
                        else:
                            self.assertTrue((expected[col] == actual[col]).all())


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
