"""Tests `dms_tools2.pacbio.CCS` class and related functions.
"""

import os
import unittest
import copy

import numpy
import pysam

import dms_tools2.pacbio


class test_pacbio_CCS(unittest.TestCase):
    """Tests `dms_tools2.pacbio.CCS` and related functions."""


    def setUp(self):
        """Initialize a `CCS` object"""

        cwd = os.path.abspath(os.path.dirname(__file__))
        indir = os.path.join(cwd, 'pacbio_ccs_files')
        reportfile = os.path.join(indir, 'test_report.txt')
        bamfile = os.path.join(indir, 'test_ccs.bam')


        self.bamlines = sum(1 for _ in pysam.AlignmentFile(bamfile,
                'rb', check_sq=False))

        self.testdir = os.path.join(cwd, 'test_pacbio_ccs_files')
        if not os.path.isdir(self.testdir):
            os.mkdir(self.testdir)

        self.ccs = dms_tools2.pacbio.CCS('test', bamfile, reportfile)


    def test_qvalsToAccuracy(self):
        """Test `qvalsToAccuracy`."""
        for i, row in self.ccs.df.iterrows():
            acc = dms_tools2.pacbio.qvalsToAccuracy(row['CCS_qvals'])
            numpy.testing.assert_allclose(acc, row['CCS_accuracy'],
                    atol=1e-5, rtol=1e-5)


    def test_zmw_report(self):
        """Test creation of `CCS.zmw_report`."""
        self.assertAlmostEqual(
                self.ccs.zmw_report.fraction.sum(),
                1,
                places=3)

        self.assertEqual(
                self.bamlines,
                (self.ccs.zmw_report
                 .query('status == "Success -- CCS generated"')
                 .number.values[0]
                 )
                )

    def test_df(self):
        """Test creation of `CCS.df`."""
        self.assertEqual(self.bamlines, self.ccs.df.shape[0])
        self.assertEqual(set(self.ccs.df.columns),
                {'name', 'CCS', 'CCS_accuracy', 'CCS_qvals',
                 'CCS_length', 'passes', 'samplename'})

    def test_plotColCorrs(self):
        """Test `plot.plotColCorrs`."""
        plotfile = os.path.join(self.testdir, 'CCSresults.pdf')
        if os.path.isfile(plotfile):
            os.remove(plotfile)
        dms_tools2.plot.plotColCorrs(self.ccs.df, plotfile,
                ['passes', 'CCS_accuracy', 'CCS_length'])
        self.assertTrue(os.path.isfile(plotfile))

    def test_summarizeCCSreports(self):
        """Test `summarizeCCSreports`."""
        plotfile = os.path.join(self.testdir, 'zmw_plot.pdf')
        if os.path.isfile(plotfile):
            os.remove(plotfile)

        ccs2 = copy.deepcopy(self.ccs)
        ccs2.samplename = 'test2'
        ccs2.df = ccs2.df.assign(samplename=ccs2.samplename)

        df = dms_tools2.pacbio.summarizeCCSreports(
                [self.ccs, ccs2], 'zmw', plotfile)

        self.assertTrue(os.path.isfile(plotfile))

        self.assertTrue(all(
                self.ccs.zmw_report.number.sort_values() ==
                df.query('sample == "test"').number.sort_values()))

        self.assertEqual(len(df), 2 * len(ccs2.zmw_report))



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
