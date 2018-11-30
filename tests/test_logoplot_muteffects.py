"""Tests ``dms2_logoplot`` with ``--muteffects`` input.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess
import random

import numpy
import pandas

import dms_tools2
import dms_tools2.prefs


class test_logoplot_muteffects(unittest.TestCase):
    """Runs ``dms2_logoplot`` with muteffects as input."""

    def setUp(self):
        """Set up input data."""

        self.testdir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './test_logoplot_muteffects_files/')
        os.makedirs(self.testdir, exist_ok=True)

        self.indir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './logoplot_input_files/')

        self.name = 'muteffects'

        # define output files
        self.outfiles = [os.path.join(self.testdir, 
                self.name + suffix) for suffix in
                ['.log', '_muteffects.pdf']]

        for f in self.outfiles:
            if os.path.isfile(f):
                os.remove(f)

    def test_logoplot_muteffects(self):
        """Run ``dms2_logoplot`` with muteffects as input."""

        prefs = pandas.read_csv(os.path.join(self.indir, 'HA_prefs.csv'))
        wts = pandas.read_csv(os.path.join(self.indir,
                              'wildtypeoverlayfile.csv'))
        muteffects = dms_tools2.prefs.prefsToMutFromWtEffects(
                        prefs,
                        dms_tools2.AAS,
                        wts)
        muteffectsfile = os.path.join(self.testdir, 'muteffects.csv')
        muteffects.to_csv(muteffectsfile, index=False)

        subprocess.check_call([
                'dms2_logoplot',
                '--name', self.name,
                '--muteffects', muteffectsfile,
                '--outdir', self.testdir,
                '--overlay1', os.path.join(self.indir, 'continuous_overlay.csv'),
                        'cp', 'continuous-property',
                '--overlay2', os.path.join(self.indir, 'discrete_overlay.csv'),
                        'dp', 'discrete-property',
                '--underlay', 'yes',
                '--scalebar', '10', '1024-fold effect',
                ])

        for f in self.outfiles:
            self.assertTrue(os.path.isfile(f), "did not create {0}".format(f))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
