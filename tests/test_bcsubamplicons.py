"""Tests ``dms2_bcsubamplicons``.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess



class test_bcsubamplicons(unittest.TestCase):
    """Runs ``dms2_bcsubamplicons`` on test data.
    """

    def setUp(self):
        """Gets files set up appropriately."""

        self.currdir = os.path.abspath(os.path.dirname(__file__))

        # defines expected files
        self.expectdir = os.path.join(self.currdir, 
                './expected_bcsubamplicons_files/')
        self.refseq = os.path.join(self.expectdir, 'WSN-HA.fasta')
        self.r1 = os.path.join(self.expectdir, 'test_R1.fastq.gz')
        self.r2 = self.r1.replace('_R1', '_R2')
        self.alignspecs = ['1,426,36,38',
                           '427,849,32,32',
                           '850,1275,31,37',
                           '1276,1698,46,45'
                          ]
        self.correctcounts = os.path.join(self.expectdir, 'correct_counts.txt')
        self.correctstats = os.path.join(self.expectdir, 
                'correct_summarystats.txt')
        for f in [self.refseq, self.r1, self.r2, 
                self.correctcounts, self.correctstats]:
            assert os.path.isfile(f), "Missing file {0}".format(f)

        # defines output files
        self.outdir = os.path.join(self.currdir,
                './test_bcsubamplicons_files/')
        self.name = 'test'
        self.counts = '{0}/{1}_counts.csv'.format(self.outdir, self.name)
        self.stats = '{0}/{1}_summarystats.csv'.format(self.outdir, self.name)
        for f in [self.counts, self.stats]:
            if os.path.isfile(f):
                os.remove(f)


    def test_dms2_bcsubamplicons(self):
        """Runs ``dms2_bcsubamplicons`` on test data."""
        cmds = [
                'dms2_bcsubamplicons',
                '--name', self.name,
                '--refseq', self.refseq,
                '--alignspecs'] + self.alignspecs + [
                '--outdir', self.outdir,
                '--R1', self.r1
               ] 
        sys.stderr.write('\nRunning the following command:\n{0}\n'.format(
                ' '.join(cmds)))
        subprocess.check_call(cmds)
#        self.assertTrue(os.path.isfile(self.counts), 
#                'Failed to create file {0}'.format(self.counts))
#        self.assertTrue(os.path.isfile(self.stats), 
#                'Failed to create file {0}'.format(self.stats))




if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
