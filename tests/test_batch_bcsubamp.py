"""Tests ``dms2_batch_bcsubamp``.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess
import random
import gzip
import numpy
import pandas



class test_batch_bcsubamp(unittest.TestCase):
    """Runs ``dms2_batch_bcsubamp`` on test data."""
    SITEMASK = []

    def setUp(self):
        """Set up input data."""

        self.testdir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './test_batch_bcsubamp_files/')
        if not os.path.isdir(self.testdir):
            os.mkdir(self.testdir)

        self.indir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './bcsubamp_input_files/')
        self.refseqfile = os.path.join(self.indir, 'WSN-HA.fasta')
        self.assertTrue(os.path.isfile(self.refseqfile))
        self.expectedcodonmuttypes = os.path.join(self.indir,
                'expected_codonmuttypes.csv')
        self.assertTrue(os.path.isfile(self.expectedcodonmuttypes))

        self.alignspecs = ['1,426,36,38', '427,849,32,32',
                           '850,1275,31,37', '1276,1698,46,45']

        self.summaryprefix = 'summary'
        if self.SITEMASK:
            self.summaryprefix += '-sitemask'

        # define output files
        files = ['.log', '_bcstats.pdf', '_codonmuttypes.pdf',
                '_codonmuttypes.csv', '_codonntchanges.pdf',
                '_depth.pdf', '_mutfreq.pdf', '_readsperbc.pdf',
                '_readstats.pdf', '_singlentchanges.pdf']
        self.outfiles = dict([(f, '{0}/{1}{2}'.format(
                self.testdir, self.summaryprefix, f)) for f in files])
        for f in self.outfiles.values():
            if os.path.isfile(f):
                os.remove(f)


    def test_dms2_batch_bcsubamp(self):
        """Runs ``dms2_batch_bcsubamp`` on test data."""

        batchfile = os.path.join(self.testdir, 'batch.csv')
        samples = ['sample-1', 'sample-2']
        df = pandas.DataFrame({'name':samples, 
                'R1':['{0}_R1.fastq.gz'.format(sample) for sample in samples]})
        df.to_csv(batchfile, index=False)

        cmds = [
                'dms2_batch_bcsubamp',
                '--batchfile', batchfile,
                '--summaryprefix', self.summaryprefix,
                '--refseq', self.refseqfile,
                '--outdir', self.testdir,
                '--fastqdir', self.indir,
                '--alignspecs',
               ] + self.alignspecs
        if self.SITEMASK:
            sitemaskfile = os.path.join(self.testdir, 'sitemask.csv')
            pandas.DataFrame({'site':self.SITEMASK, 
                    'extracol':['a'] * len(self.SITEMASK)}).to_csv(sitemaskfile)
            cmds += ['--sitemask', sitemaskfile]
        sys.stderr.write('\nRunning the following command:\n{0}\n'.format(
                ' '.join(cmds)))
        subprocess.check_call(cmds)

        for f in self.outfiles.values():
            self.assertTrue(os.path.isfile(f), "Failed to create {0}".format(f))

        if self.SITEMASK:
            # just test right site columns 
            for sample in samples:
                counts = pandas.read_csv(os.path.join(self.testdir,
                        sample + '_codoncounts.csv'))
                self.assertTrue(all(counts['site'] == self.SITEMASK))

        else:
            # test for correct output
            expected = pandas.read_csv(self.expectedcodonmuttypes, 
                    index_col='name')
            actual = pandas.read_csv(self.outfiles['_codonmuttypes.csv'], 
                    index_col='name')
            for c in expected.columns:
                self.assertTrue(numpy.allclose(expected[c].values, actual[c].values))


class test_batch_bcsubamp_sitemask(test_batch_bcsubamp):
    """Runs ``dms2_batch_bcsubamp`` with ``--sitemask`` on test data."""
    SITEMASK = list(range(2, 50))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
