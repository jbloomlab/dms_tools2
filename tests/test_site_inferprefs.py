"""Tests inference of preferences for a site."""


import sys
import os
import unittest
import random
import numpy
import dms_tools2.prefs
from dms_tools2 import AAS, NTS, CODONS



class TestInferSitePreferences(unittest.TestCase):
    """Tests preference inference when priors are exactly correct.

    Makes sure inference converges to values that are very close
    to the correct one for all three error models:

        1) No errors

        2) Same errors pre- and post-selection

        3) Different errors pre- and post-selection

    Then also makes sure that use a no-error model for data with
    errors does **not** converge.
    """
    SEED = 1
    CHARLIST = NTS

    def setUp(self):
        """Configuration information for test."""
        self.n_jobs = -1

    def test_inferSitePrefs(self):
        """Inference with with correct priors."""
        # rates drawn uniformly from following ranges after dividing by nchars
        mutrate = (0.005, 0.025) 
        errorratepre = (0.0001, 0.0003) 
        errorratepost = (0.0002, 0.0005) 

        # prefs drawn uniformly from this range and adjusted to sum to one
        pidist = (1e-5, 0.6) 

        random.seed(self.SEED)
        numpy.random.seed(self.SEED)
        wtchar = random.choice(self.CHARLIST)
        iwtchar = self.CHARLIST.index(wtchar)
        nchars = len(self.CHARLIST)
        mur = []
        pir = []
        epsilonr = []
        rhor = []
        for i in range(nchars):
            mur.append(random.uniform(mutrate[0] / float(nchars), 
                    mutrate[1] / float(nchars)))
            epsilonr.append(random.uniform(errorratepre[0] / float(nchars), 
                    errorratepre[1] / float(nchars)))
            rhor.append(random.uniform(errorratepost[0] / float(nchars), 
                    errorratepost[1] / float(nchars)))
            pir.append(random.uniform(pidist[0], pidist[1]))
        mur[iwtchar] = 1.0 - sum([mur[i] for i in range(nchars)
                if i != iwtchar])
        epsilonr[iwtchar] = 1.0 - sum([epsilonr[i] for i in range(nchars) 
                if i != iwtchar])
        rhor[iwtchar] = 1.0 - sum([rhor[i] for i in range(nchars) 
                if i != iwtchar])
        pir[iwtchar] = 1.0
        pir = numpy.array(pir)
        pir /= pir.sum()
        mur = numpy.array(mur)
        epsilonr = numpy.array(epsilonr)
        rhor = numpy.array(rhor)
        deltar = numpy.zeros(nchars)
        deltar[iwtchar] = 1.0
        priors = {
                'pir_prior_params':dict([(char, pir[i] * nchars) 
                        for (i, char) in enumerate(self.CHARLIST)]),
                'mur_prior_params':dict([(char, mur[i] * nchars) 
                        for (i, char) in enumerate(self.CHARLIST)]),
                'epsilonr_prior_params':dict([(char, epsilonr[i] * nchars) 
                        for (i, char) in enumerate(self.CHARLIST)]),
                'rhor_prior_params':dict([(char, rhor[i] * nchars) 
                        for (i, char) in enumerate(self.CHARLIST)]),
                }
        difflist = []
        for (depth, maxdiffsum) in [(1e9, 0.01)]:

            # first with no errors
            nrpre = numpy.random.multinomial(depth, mur)
            nrpost = numpy.random.multinomial(depth, mur * pir / 
                    numpy.dot(mur, pir))
            counts = {
                    'pre':dict([(char, nrpre[i]) for (i, char) 
                            in enumerate(self.CHARLIST)]),
                    'post':dict([(char, nrpost[i]) for (i, char) 
                            in enumerate(self.CHARLIST)]),
            }
            (converged, pi_means, pi_95credint, logstring) = \
                    dms_tools2.prefs.inferSitePrefs(self.CHARLIST, wtchar, 
                    'none', counts, priors, n_jobs=self.n_jobs)
            self.assertTrue(converged, 'failed to converge')
            inferred_pi = [pi_means[char] for char in self.CHARLIST]
            diffsum = sum([abs(x - y) for (x, y) in zip(pir, inferred_pi)])
            self.assertTrue(diffsum < maxdiffsum, 'wrong preferences')

            # now with same errors pre and post
            nrpre = numpy.random.multinomial(depth, mur + epsilonr - deltar)
            nrpost = numpy.random.multinomial(depth, mur * pir / 
                    numpy.dot(mur, pir) + epsilonr - deltar)
            nrerr = numpy.random.multinomial(depth, epsilonr)
            counts = {
                    'pre':dict([(char, nrpre[i]) for (i, char) 
                            in enumerate(self.CHARLIST)]),
                    'post':dict([(char, nrpost[i]) for (i, char) 
                            in enumerate(self.CHARLIST)]),
                    'err':dict([(char, nrerr[i]) for (i, char) 
                            in enumerate(self.CHARLIST)]),
                    }
            (converged, pi_means, pi_95credint, logstring) = \
                    dms_tools2.prefs.inferSitePrefs(self.CHARLIST, wtchar, 
                    'same', counts, priors, n_jobs=self.n_jobs)
            self.assertTrue(converged, 'failed to converge')
            inferred_pi = [pi_means[char] for char in self.CHARLIST]
            diffsum = sum([abs(x - y) for (x, y) in zip(pir, inferred_pi)])
            self.assertTrue(diffsum < maxdiffsum, 'wrong preferences')

            # now with different errors pre and post
            nrpre = numpy.random.multinomial(depth, mur + epsilonr - deltar)
            nrpost = numpy.random.multinomial(depth, mur * pir / 
                    numpy.dot(mur, pir) + rhor - deltar)
            nrerrpre = numpy.random.multinomial(depth, epsilonr)
            nrerrpost = numpy.random.multinomial(depth, rhor)
            counts = {
                    'pre':dict([(char, nrpre[i]) for (i, char) 
                            in enumerate(self.CHARLIST)]),
                    'post':dict([(char, nrpost[i]) for (i, char) 
                            in enumerate(self.CHARLIST)]),
                    'errpre':dict([(char, nrerrpre[i]) for (i, char) 
                            in enumerate(self.CHARLIST)]),
                    'errpost':dict([(char, nrerrpost[i]) for (i, char) 
                            in enumerate(self.CHARLIST)]),
                    }
            (converged, pi_means, pi_95credint, logstring) = \
                    dms_tools2.prefs.inferSitePrefs(self.CHARLIST, wtchar, 
                    'different', counts, priors, n_jobs=self.n_jobs)
            self.assertTrue(converged, 'failed to converge')
            inferred_pi = [pi_means[char] for char in self.CHARLIST]
            diffsum = sum([abs(x - y) for (x, y) in zip(pir, inferred_pi)])
            self.assertTrue(diffsum < maxdiffsum, 'wrong preferences')

            # make sure does NOT converge without error estimates
            (converged, pi_means, pi_95credint, logstring) = \
                    dms_tools2.prefs.inferSitePrefs(self.CHARLIST, wtchar, 
                    'none', counts, priors, n_jobs=self.n_jobs)
            self.assertTrue(converged, 'prefs failed to converge')
            inferred_pi = [pi_means[char] for char in self.CHARLIST]
            diffsum = sum([abs(x - y) for (x, y) in zip(pir, inferred_pi)])
            self.assertFalse(diffsum < maxdiffsum, 'incorrectly converged')


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
