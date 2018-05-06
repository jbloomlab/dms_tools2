"""Tests `dms_tools2.pacbio.CCS` class and related functions.

Requires ``samtools`` to be installed.
"""

import os
import unittest
import copy

import numpy
import pandas
from pandas.testing import assert_frame_equal

import dms_tools2.pacbio


class test_pacbio_CCS(unittest.TestCase):
    """Tests `dms_tools2.pacbio.CCS` and related functions."""


    def setUp(self):
        """Initialize a `CCS` object"""

        cwd = os.path.abspath(os.path.dirname(__file__))
        indir = os.path.join(cwd, 'pacbio_ccs_files')
        reportfile = os.path.join(indir, 'test_report.txt')
        bamfile = os.path.join(indir, 'test_ccs.bam')

        self.testdir = os.path.join(cwd, 'test_pacbio_ccs_files')
        if not os.path.isdir(self.testdir):
            os.mkdir(self.testdir)

        self.ccs = dms_tools2.pacbio.CCS('test', bamfile, reportfile)


    def test_samfile(self):
        """Test creation of `CCS.samfile`."""
        samfile = self.ccs.samfile
        self.assertTrue(os.path.isfile(samfile))


    def test_zmw_report(self):
        """Test creation of `CCS.zmw_report`."""
        self.assertAlmostEqual(
                self.ccs.zmw_report.fraction.sum(),
                1,
                places=3)

        with open(self.ccs.samfile) as f:
            samlines = len(f.readlines())
        self.assertEqual(
                samlines,
                (self.ccs.zmw_report
                 .query('status == "Success -- CCS generated"')
                 .number.values[0]
                 )
                )

    def test_summarizeCCSreports(self):
        """Test `summarizeCCSreports`."""
        plotfile = os.path.join(self.testdir, 'zmw_plot.pdf')
        if os.path.isfile(plotfile):
            os.remove(plotfile)

        ccs2 = copy.deepcopy(self.ccs)
        ccs2.name = 'test2'

        df = dms_tools2.pacbio.summarizeCCSreports(
                [self.ccs, ccs2], 'zmw', plotfile)

        self.assertTrue(os.path.isfile(plotfile))

        self.assertTrue(all(
                self.ccs.zmw_report.number.sort_values() ==
                df.query('sample == "test"').number.sort_values()))

        self.assertEqual(len(df), 2 * len(ccs2.zmw_report))


    def tearDown(self):
        """Remove created files."""
        try:
            os.remove(self.ccs.samfile)
        except:
            pass



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
