"""Tests `dms_tools2.prefs.inferPrefsByRatio`."""


import sys
import os
import unittest
import random
import numpy
import pandas
import dms_tools2
import dms_tools2.prefs



class TestInferPrefsByRatio(unittest.TestCase):
    """Tests `dms_tools2.inferPrefsByRatio`.
    """

    #: draw Nr from this range for different samples
    Nr_RANGE = (4000, 4000)
    # tolerance for comparisons is stringent for this case
    ATOL = 1e-8
    RTOL = 1e-5

    def setUp(self):
        """Set up for tests."""
        random.seed(1)

        self.charlist = dms_tools2.AAS
        self.nchars = len(self.charlist)
        self.sites = ['1', '1A', '2']
        self.nsites = len(self.sites)
        self.wts = [random.choice(self.charlist) for r in self.sites]
        self.pseudocount = 1

    def test_inferPrefsByRatio_NoCounts(self):
        """Prefs should be uniform with no counts."""
        dfs = {}
        for ctype in ['pre', 'post', 'errpre', 'errpost']:
            dfs[ctype] = pandas.DataFrame(dict([('site', self.sites)]
                    + [(c, [0] * self.nsites) for c in self.charlist]))
        for (errpre, errpost) in [
                (dfs['errpre'], dfs['errpost']),
                (None, None),
                ]:
            prefs = dms_tools2.prefs.inferPrefsByRatio(self.charlist, 
                    self.sites, self.wts, dfs['pre'], dfs['post'], 
                    errpre, errpost, self.pseudocount)
            for c in self.charlist:
                self.assertTrue(numpy.allclose(prefs[c].values,
                        numpy.full(self.nsites, 1.0 / self.nchars)),
                        str(prefs))

    def test_inferPrefsByRatio_Simplecounts(self):
        """Some simple non-zero counts."""

        # simulate some counts
        random.seed(1)
        numpy.random.seed(1)
        dfs = {}
        for ctype in ['pre', 'post', 'errpre', 'errpost']:
            rows = []
            for (r, wt) in zip(self.sites, self.wts):
                Nr = random.uniform(*self.Nr_RANGE)
                pconc = numpy.ones(self.nchars)
                if 'err' in ctype:
                    pconc[self.charlist.index(wt)] = 25
                else:
                    pconc[self.charlist.index(wt)] = 5
                p = numpy.random.dirichlet(pconc)
                counts = numpy.random.multinomial(Nr, p)
                rows.append([r] + list(counts))
            dfs[ctype] = pandas.DataFrame(rows, 
                    columns=['site'] + self.charlist)

        # compute prefs using function
        prefs_noerr = dms_tools2.prefs.inferPrefsByRatio(self.charlist,
                self.sites, self.wts, dfs['pre'], dfs['post'],
                None, None, self.pseudocount)
        self.assertTrue(list(prefs_noerr['site']) == self.sites)
        prefs_err = dms_tools2.prefs.inferPrefsByRatio(self.charlist,
                self.sites, self.wts, dfs['pre'], dfs['post'],
                dfs['errpre'], dfs['errpost'], self.pseudocount)
        self.assertTrue(list(prefs_err['site']) == self.sites)

        # compare to values calculated by hand following notation
        # in docs for `dms_tools2.prefs.inferPrefsByRatio`
        P = float(self.pseudocount) 
        A = self.nchars
        for (r, wt) in zip(self.sites, self.wts):
            f_before_noerr = numpy.ndarray(self.nchars)
            f_after_noerr = numpy.ndarray(self.nchars)
            f_before_err = numpy.ndarray(self.nchars)
            f_after_err = numpy.ndarray(self.nchars)
            for (i, a) in enumerate(self.charlist):
                if a == wt:
                    iwt = i
                    delta = 1
                else:
                    delta = 0
                f = {}
                Nr = {}
                for ctype in ['pre', 'post', 'errpre', 'errpost']:
                    df = dfs[ctype]
                    nr = df[df['site'] == r][self.charlist]
                    Nr[ctype] = nr.values.sum()
                    nra = nr[a].values[0]
                    f[ctype] = (nra + P) / (Nr[ctype] + A * P)
                f_before_noerr[i] = max(P / (Nr['pre'] + A * P), f['pre'])
                f_after_noerr[i] = max(P / (Nr['post'] + A * P), f['post'])
                f_before_err[i] = max(P / (Nr['pre'] + A * P), 
                        f['pre'] + delta - f['errpre'])
                f_after_err[i] = max(P / (Nr['post'] + A * P), 
                        f['post'] + delta - f['errpost'])
            phi_noerr = (f_after_noerr / f_after_noerr[iwt]) / (
                    f_before_noerr / f_before_noerr[iwt])
            phi_err = (f_after_err / f_after_err[iwt]) / (
                    f_before_err / f_before_err[iwt])
            pi_noerr = phi_noerr / phi_noerr.sum()
            pi_err = phi_err / phi_err.sum()
            self.assertTrue(numpy.allclose(pi_noerr,
                    prefs_noerr.query('site == @r')[self.charlist],
                    atol=self.ATOL, rtol=self.RTOL))
            self.assertFalse(numpy.allclose(pi_err,
                    prefs_noerr.query('site == @r')[self.charlist],
                    atol=self.ATOL, rtol=self.RTOL))
            self.assertTrue(numpy.allclose(pi_err,
                    prefs_err.query('site == @r')[self.charlist],
                    atol=self.ATOL, rtol=self.RTOL))
            self.assertFalse(numpy.allclose(pi_noerr,
                    prefs_err.query('site == @r')[self.charlist],
                    atol=self.ATOL, rtol=self.RTOL))



class TestInferPrefsByRatioDiffDepth(TestInferPrefsByRatio):
    """Tests `dms_tools2.inferPrefsByRatio` with different depths.

    There is pseudocount scaling in this case. We handle this by
    calculating expected results without scaling and just relaxing
    tolerance while keeping depth close among samples.
    """

    #: draw Nr from this range for different samples
    Nr_RANGE = (300000, 400000)
    # tolerance for comparisons is stringent for this case
    ATOL = 1e-3
    RTOL = 1e-3


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
