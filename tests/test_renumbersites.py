"""Tests `dms_tools2.utils.renumberSites`."""


import os
import unittest
import numpy
import pandas
from dms_tools2.utils import renumberSites


class test_renumberSites(unittest.TestCase):
    """Tests `dms_tools2.utils.renumberSites`."""

    def test_renumberSites(self):
        """Tests `dms_tools2.utils.renumberSites`."""

        testdir = os.path.join(os.path.dirname(__file__),
                'test_renumbersites_files')
        if not os.path.isdir(testdir):
            os.mkdir(testdir)

        renumbfile = os.path.join(testdir, 'renumbfile.csv')
        with open(renumbfile, 'w') as f:
            f.write('\n'.join([
                    'original,new',
                    '1,3',
                    '2,3a',
                    '3,None',
                    '4,4',
                    ]))

        infile1 = os.path.join(testdir, 'file1.csv')
        infile2 = os.path.join(testdir, 'file2.csv')
        with open(infile1, 'w') as f:
            f.write('\n'.join([
                    'site,val',
                    '1,A',
                    '2,B',
                    '3,C',
                    '4,D',
                    '5,E',
                    ]))
        with open(infile2, 'w') as f:
            f.write('\n'.join([
                    'site,val',
                    '1,A',
                    '2,B',
                    '2,B',
                    ]))

        outdir = os.path.join(testdir, 'outdir')
        outfiles = [os.path.join(outdir, 'file1.csv'),
                    os.path.join(outdir, 'file2.csv')]
        for f in outfiles:
            if os.path.isfile(f):
                os.remove(f)

        with self.assertRaises(ValueError):
            renumberSites(renumbfile, [infile1, infile2],
                    missing='error', outdir=outdir)

        renumberSites(renumbfile, [infile1, infile2],
                missing='skip', outdir=outdir)
        self.assertTrue(all(map(os.path.isfile, outfiles)))
        df1 = pandas.read_csv(outfiles[0])
        self.assertTrue(all(df1['site'] == ['3', '3a', '4', '5']))
        self.assertTrue(all(df1['val'] == ['A', 'B', 'D', 'E']))
        df2 = pandas.read_csv(outfiles[1])
        self.assertTrue(all(df2['site'] == ['3', '3a', '3a']))
        self.assertTrue(all(df2['val'] == ['A', 'B', 'B']))

        for f in outfiles:
            os.remove(f)

        renumberSites(renumbfile, [infile1, infile2],
                missing='drop', outfiles=outfiles)
        self.assertTrue(all(map(os.path.isfile, outfiles)))
        df1 = pandas.read_csv(outfiles[0])
        self.assertTrue(all(df1['site'] == ['3', '3a', '4']))
        self.assertTrue(all(df1['val'] == ['A', 'B', 'D']))
        df2 = pandas.read_csv(outfiles[1])
        self.assertTrue(all(df2['site'] == ['3', '3a', '3a']))
        self.assertTrue(all(df2['val'] == ['A', 'B', 'B']))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
