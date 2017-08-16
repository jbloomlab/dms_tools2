"""
========================
prefs
========================

Performs operations related to estimating site-specific amino-acid
preferences.

Uses `pystan <https://pystan.readthedocs.io/en/latest>`_
to perform MCMC for Bayesian inferences.
"""


import sys
import tempfile
import time
import math
import pickle
import numpy
import numpy.random
import pystan

#: minimum value for Dirichlet prior elements
PRIOR_MIN_VALUE = 1.0e-7 


class StanModelNoneErr(object):
    """``pystan`` model when `error_model` is `none`.
    
    For use by inferSitePrefs`."""
    def __init__(self):
        self.pystancode =\
"""
data {
    int<lower=1> Nchar; // 64 for codons, 20 for amino acids, 4 for nucleotides
    int<lower=0> nrpre[Nchar]; // counts pre-selection
    int<lower=0> nrpost[Nchar]; // counts post-selection
    vector<lower=%g>[Nchar] pir_prior_params; // Dirichlet prior params
    vector<lower=%g>[Nchar] mur_prior_params; // Dirichlet prior params
}
parameters {
    simplex[Nchar] pir;
    simplex[Nchar] mur;
}
transformed parameters {
    simplex[Nchar] fr;
    fr <- pir .* mur / dot_product(pir, mur);
}
model {
    pir ~ dirichlet(pir_prior_params);
    mur ~ dirichlet(mur_prior_params);
    nrpre ~ multinomial(mur);
    nrpost ~ multinomial(fr);
}
""" % (PRIOR_MIN_VALUE, PRIOR_MIN_VALUE)
        self.model = pystan.StanModel(model_code=self.pystancode)


class StanModelSameErr:
    """``pystan`` model when `error_model` is `same`.
    
    For use by inferSitePrefs`."""
    def __init__(self):
        self.pystancode =\
"""
data {
    int<lower=1> Nchar; // 64 for codons, 20 for amino acids, 4 for nucleotides
    int<lower=1, upper=Nchar> iwtchar; // index of wildtype character in 1, ... numbering
    int<lower=0> nrpre[Nchar]; // counts pre-selection
    int<lower=0> nrpost[Nchar]; // counts post-selection
    int<lower=0> nrerr[Nchar]; // counts in error control
    vector<lower=%g>[Nchar] pir_prior_params; // Dirichlet prior params
    vector<lower=%g>[Nchar] mur_prior_params; // Dirichlet prior params
    vector<lower=%g>[Nchar] epsilonr_prior_params; // Dirichlet prior params
}
transformed data {
    simplex[Nchar] deltar;
    for (ichar in 1:Nchar) {
        deltar[ichar] <- 0.0;
    }
    deltar[iwtchar] <- 1.0;
}
parameters {
    simplex[Nchar] pir;
    simplex[Nchar] mur;
    simplex[Nchar] epsilonr;
}
transformed parameters {
    simplex[Nchar] fr_plus_err;
    simplex[Nchar] mur_plus_err;
    fr_plus_err <- pir .* mur / dot_product(pir, mur) + epsilonr - deltar;
    mur_plus_err <- mur + epsilonr - deltar;
}
model {
    pir ~ dirichlet(pir_prior_params);
    mur ~ dirichlet(mur_prior_params);
    epsilonr ~ dirichlet(epsilonr_prior_params);
    nrerr ~ multinomial(epsilonr);
    nrpre ~ multinomial(mur_plus_err);
    nrpost ~ multinomial(fr_plus_err);
}
""" % (PRIOR_MIN_VALUE, PRIOR_MIN_VALUE, PRIOR_MIN_VALUE)
        self.model = pystan.StanModel(model_code=self.pystancode)


class StanModelDifferentErr:
    """``pystan`` model when `error_model` is `different`.
    
    For use by inferSitePrefs`."""
    def __init__(self):
        self.pystancode =\
"""
data {
    int<lower=1> Nchar; // 64 for codons, 20 for amino acids, 4 for nucleotides
    int<lower=1, upper=Nchar> iwtchar; // index of wildtype character in 1, ... numbering
    int<lower=0> nrpre[Nchar]; // counts pre-selection
    int<lower=0> nrpost[Nchar]; // counts post-selection
    int<lower=0> nrerrpre[Nchar]; // counts in pre-selection error control
    int<lower=0> nrerrpost[Nchar]; // counts in post-selection error control
    vector<lower=%g>[Nchar] pir_prior_params; // Dirichlet prior params
    vector<lower=%g>[Nchar] mur_prior_params; // Dirichlet prior params
    vector<lower=%g>[Nchar] epsilonr_prior_params; // Dirichlet prior params
    vector<lower=%g>[Nchar] rhor_prior_params; // Dirichlet prior params
}
transformed data {
    simplex[Nchar] deltar;
    for (ichar in 1:Nchar) {
        deltar[ichar] <- 0.0;
    }
    deltar[iwtchar] <- 1.0;
}
parameters {
    simplex[Nchar] pir;
    simplex[Nchar] mur;
    simplex[Nchar] epsilonr;
    simplex[Nchar] rhor;
}
transformed parameters {
    simplex[Nchar] fr_plus_err;
    simplex[Nchar] mur_plus_err;
    fr_plus_err <- pir .* mur / dot_product(pir, mur) + rhor - deltar;
    mur_plus_err <- mur + epsilonr - deltar;
}
model {
    pir ~ dirichlet(pir_prior_params);
    mur ~ dirichlet(mur_prior_params);
    epsilonr ~ dirichlet(epsilonr_prior_params);
    rhor ~ dirichlet(rhor_prior_params);
    nrerrpre ~ multinomial(epsilonr);
    nrerrpost ~ multinomial(rhor);
    nrpre ~ multinomial(mur_plus_err);
    nrpost ~ multinomial(fr_plus_err);
}
""" % (PRIOR_MIN_VALUE, PRIOR_MIN_VALUE, PRIOR_MIN_VALUE, PRIOR_MIN_VALUE)
        self.model = pystan.StanModel(model_code=self.pystancode)


def _InitialValuePreferences(error_model, nchains, iwtchar, nchars):
    """Gets valid initial values for PyStan preference inference.

    The values automatically initialized by PyStan frequently have invalid values.
    This function will generate random values that are valid
    for initialization and return a list that can be passed to the *StanModel*
    as the *init* argument.
    """
    initializationattempts = 10 # multiple attempts, as might fail for pathological random values
    nrescales = 10 # rescale non-wildtype values down this many times
    rescalefactor = 5.0 # rescale non-wildtype values down by this much
    deltar = numpy.zeros(nchars)
    deltar[iwtchar] = 1.0
    init = []
    concentrationparameters = [1.0 for ichar in range(nchars)]
    concentrationparameters[iwtchar] = 10.0 # make this element bigger, as it probably should be
    if error_model == 'none':
        for chain in range(nchains):
            while True:
                chain_init = {\
                    'pir':numpy.random.dirichlet(concentrationparameters),\
                    'mur':numpy.random.dirichlet(concentrationparameters),\
                    }
                if all(chain_init['pir'] > PRIOR_MIN_VALUE) and all(chain_init['mur'] > PRIOR_MIN_VALUE):
                    init.append(chain_init)
                    break
        return init
    elif error_model == 'same':
        posterrname = 'epsilonr'
    elif error_model == 'different':
        posterrname = 'rhor'
    else:
        raise ValueError("Invalid error_model of %s when getting initial values" % error_model)
    for chain in range(nchains):
        for iattempt in range(initializationattempts):
            chain_init = {\
                    'pir':numpy.random.dirichlet(concentrationparameters),\
                    'mur':numpy.random.dirichlet(concentrationparameters),\
                    'epsilonr':numpy.random.dirichlet(concentrationparameters),\
                    }
            if error_model == 'different':
                chain_init['rhor'] = numpy.random.dirichlet(concentrationparameters)
            irescale = 0
            while irescale < nrescales and (any(chain_init['pir'] * chain_init['mur'] / numpy.dot(chain_init['pir'], chain_init['mur']) + chain_init[posterrname] - deltar < PRIOR_MIN_VALUE) or any(chain_init['mur'] + chain_init['epsilonr'] - deltar < PRIOR_MIN_VALUE)): 
                irescale += 1
                chain_init['epsilonr'] /= rescalefactor
                chain_init['epsilonr'][iwtchar] = 1.0 - sum([chain_init['epsilonr'][ichar] for ichar in range(nchars) if ichar != iwtchar])
                if error_model == 'different':
                    chain_init['rhor'] /= rescalefactor
                    chain_init['rhor'][iwtchar] = 1.0 - sum([chain_init['rhor'][ichar] for ichar in range(nchars) if ichar != iwtchar])
            if all(chain_init['pir'] * chain_init['mur'] / numpy.dot(chain_init['pir'], chain_init['mur']) + chain_init[posterrname] - deltar > PRIOR_MIN_VALUE) and all(chain_init['mur'] + chain_init['epsilonr'] - deltar > PRIOR_MIN_VALUE):
                break
        else:
            raise ValueError("Failed to initialize for same even after %d attempts" % initializationattempts)
        init.append(chain_init)
    return init


def InferSitePreferencesFromEnrichmentRatios(characterlist, wtchar, error_model, counts, pseudocounts=1):
    r"""Infers site-specific preferences from enrichment ratios.

    This function mirrors the operations performed by *InferSitePreferences*, expect the preferences
    are calculated directly from enrichment ratios. 

    *characterlist*, *wtchar*, *error_model*, and *counts* have the same meaning as for *InferSitePreferences*.

    *pseudocounts* is a number > 0. If the counts for a character are less than *pseudocounts*, either due to low counts or error correction, then the counts for that character are changed to *pseudocounts* to avoid estimating ratios of zero, less than zero, or infinity.

    Briefly, we set the enrichment ratio of the wildtype character :math:`\rm{wt}` at this site equal to one
    
    .. math::
       
        \phi_{wt} = 1
    
    Then, for each non-wildtype character :math:`x`, we calculate the enrichment relative to :math:`\rm{wt}` as

    .. math::
       
        \phi_{x} = 
        \frac{
          \left(\frac{
            \max\left(\frac{\mathcal{P}}{N_r^{\rm{post}}}, \frac{n_{r,x}^{\rm{post}}}{N_r^{\rm{post}}} - \frac{n_{r,x}^{\rm{errpost}}}{N_r^{\rm{errpost}}}\right)
          }{
            \frac{n_{r,\rm{wt}}^{\rm{post}}}{N_r^{\rm{post}}} + \delta - \frac{n_{r,\rm{wt}}^{\rm{errpost}}}{N_r^{\rm{errpost}}}
          }\right)
        }
        {
          \left(\frac{
            \max\left(\frac{\mathcal{P}}{N_r^{\rm{pre}}}, \frac{n_{r,x}^{\rm{pre}}}{N_r^{\rm{pre}}} - \frac{n_{r,x}^{\rm{errpre}}}{N_r^{\rm{errpre}}}\right)
          }{
            \frac{n_{r,\rm{wt}}^{\rm{pre}}}{N_r^{\rm{pre}}} + \delta - \frac{n_{r,\rm{wt}}^{\rm{errpre}}}{N_r^{\rm{errpre}}}
          }\right)
        }

    where :math:`\mathcal{P}` is the value of *pseudocounts*. When *error_model* is *none*, then all terms involving the error corrections (with superscript *err*) are ignored and :math:`\delta` is set to zero; otherwise :math:`\delta` is one.

    Next, for each character :math:`x`, including :math:`\rm{wt}`, we calculate the preference for :math:`x` as

    .. math::

        \pi_x = \frac{\phi_x}{\sum_y \phi_y}.

    The return value is: *(converged, pi_means, pi_95credint, logstring)*, where the tuple
    entries have the same meaning as for *InferSitePreferences* except that *pi_95credint* is
    *None* since no credible intervals can be estimated from direct enrichment ratio calculation
    as it is not a statistical model, and *converged* is *True* since this calculation
    always converges.

    For testing the code using doctest:
    
    >>> # Hypothetical data for a site
    >>> characterlist = ['A', 'T', 'G', 'C']
    >>> wtchar = 'A'
    >>> counts = {}
    >>> counts['nrpost'] = {'A':310923, 'T':13, 'C':0, 'G':37}
    >>> counts['nrerrpost'] = {'A':310818, 'T':0, 'C':0, 'G':40}
    >>> counts['nrerr'] = {'A':310818, 'T':0, 'C':0, 'G':40} # Same as 'nrerrpost'
    >>> counts['nrpre'] = {'A':390818, 'T':50, 'C':0, 'G':80}
    >>> counts['nrerrpre'] = {'A':390292, 'T':0, 'C':5, 'G':9}

    >>> # Using error_model = 'none'
    >>> (converged, pi, pi95, logstring) = InferSitePreferencesFromEnrichmentRatios(characterlist, wtchar, 'none', counts, pseudocounts=1)
    >>> [round(pi[x], 9) for x in characterlist]
    [0.315944301, 0.10325369, 0.18367243, 0.397129578]

    >>> # Using error_model = 'same'
    >>> (converged, pi, pi95, logstring) = InferSitePreferencesFromEnrichmentRatios(characterlist, wtchar, 'same', counts, pseudocounts=1)
    >>> [round(pi[x], 9) for x in characterlist]
    [0.380792732, 0.124446795, 0.016118953, 0.47864152]

    >>> # Using error_model = 'different'
    >>> (converged, pi, pi95, logstring) = InferSitePreferencesFromEnrichmentRatios(characterlist, wtchar, 'different', counts, pseudocounts=1)
    >>> [round(pi[x], 9) for x in characterlist]
    [0.384418849, 0.125620184, 0.006806413, 0.483154554]
    """
    
    assert pseudocounts > 0, "pseudocounts must be greater than zero, invalid value of %g" % pseudocounts
    assert wtchar in characterlist, "wtchar %s not in characterlist %s" % (wtchar, str(characterlist))
    logstring = '\tComputed preferences directly from enrichment ratios.'

    Nrpost = 0.0
    Nrpre = 0.0
    Nrerrpost = 0.0
    Nrerrpre = 0.0
    for y in characterlist:
        Nrpost += counts['nrpost'][y]
        Nrpre += counts['nrpre'][y]     
        if error_model == 'same':
            Nrerrpost += counts['nrerr'][y]
            Nrerrpre += counts['nrerr'][y]
        elif error_model == 'different':
            Nrerrpost += counts['nrerrpost'][y]
            Nrerrpre += counts['nrerrpre'][y]
        elif error_model != 'none':
            raise ValueError("Invalid error_model of %s" % error_model)
    
    psi = {wtchar:1.0}
    nrpostwt = counts['nrpost'][wtchar]
    nrprewt = counts['nrpre'][wtchar]
    for x in characterlist:
        if x == wtchar:
            continue
        nrpost = counts['nrpost'][x]
        nrpre = counts['nrpre'][x]
        if error_model == 'none':
            nrerrpre = nrerrpost = 0.0
            nrerrprewt = nrerrpostwt = 0.0
            Nrerrpre = Nrerrpost = 1.0 # Set equal to pseudocount of 1.0 to avoid dividing by zero; however, the psuedocount will make no difference in the end since all fractions with these variables will end up equalling zero anyways.
            delta = 0.0
        elif error_model == 'same':
            nrerrpre = nrerrpost = counts['nrerr'][x]
            nrerrprewt = nrerrpostwt = counts['nrerr'][wtchar]
            assert Nrerrpost == Nrerrpre
            delta = 1.0
        elif error_model == 'different':
            nrerrpre = counts['nrerrpre'][x]
            nrerrpost = counts['nrerrpost'][x]
            nrerrprewt = counts['nrerrpre'][wtchar]
            nrerrpostwt = counts['nrerrpost'][wtchar]
            delta = 1.0
        else:
            raise ValueError("Invalid error_model of %s" % error_model)
        postratio = max(pseudocounts/Nrpost, (nrpost/Nrpost)-(nrerrpost/Nrerrpost))/((nrpostwt/Nrpost)+delta-(nrerrpostwt/Nrerrpost))
        preratio = max(pseudocounts/Nrpre, (nrpre/Nrpre)-(nrerrpre/Nrerrpre))/((nrprewt/Nrpre)+delta-(nrerrprewt/Nrerrpre))
        psi[x] = postratio / preratio
    
    assert abs(psi[wtchar] - 1) < 1.0e-5, "wtchar does not have enrichment ratio of one: %g" % (psi[wtchar])
    denom = sum(psi.values())
    pi = dict([(x, psi[x] / float(denom)) for x in characterlist])
    return (True, pi, None, logstring)
    

def inferSitePrefs(characterlist, wtchar, error_model, counts, 
        priors, seed=1, niter=10000, increasetries=6, n_jobs=1, 
        r_max=1.1, neff_min=100, nchains=4, increasefactor=2):
    """Infers site-specific preferences by MCMC for a specific site.

    Infer the site-specific preferences :math:`\pi_{r,a}` for some site
    :math:`r` for each character :math:`a` by integrating over the 
    posterior defined by the product of the following priors and 
    likelihoods, where :math:`\\boldsymbol{\mathbf{\pi_r}}` indicates
    the vector of :math:`\pi_{r,a}` values for all characters:

    .. math::

        \Pr\left(\\boldsymbol{\mathbf{\pi_r}}\\right) 
        &= \mathrm{Dirichlet}\left(\\boldsymbol{\mathbf{a_{\pi,r}}}\\right)

        \Pr\left(\\boldsymbol{\mathbf{\mu_r}}\\right) 
        &= \mathrm{Dirichlet}\left(\\boldsymbol{\mathbf{a_{\mu,r}}}\\right)

        \Pr\left(\\boldsymbol{\mathbf{\epsilon_r}}\\right) 
        &= \mathrm{Dirichlet}\left(\\boldsymbol{\mathbf{a_{\epsilon,r}}}\\right)

        \Pr\left(\\boldsymbol{\mathbf{\\rho_r}}\\right) 
        &= \mathrm{Dirichlet}\left(\\boldsymbol{\mathbf{a_{\\rho,r}}}\\right)

        \Pr\left(\\boldsymbol{\mathbf{n_r^{\\rm{pre}}}} \mid \\boldsymbol{\mathbf{\mu_r}}, \\boldsymbol{\mathbf{\epsilon_r}}\\right) 
        &= \mathrm{Multinomial}\left(\\boldsymbol{\mathbf{n_r^{\\rm{pre}}}}; \\boldsymbol{\mathbf{\mu_r}} + \\boldsymbol{\mathbf{\epsilon_r}} - \\boldsymbol{\mathbf{\delta_r}}\\right)

        \Pr\left(\\boldsymbol{\mathbf{n_r^{\\rm{post}}}} \mid \\boldsymbol{\mathbf{\mu_r}}, \\boldsymbol{\mathbf{\epsilon_r}}\\right) 
        &= \mathrm{Multinomial}\left(\\boldsymbol{\mathbf{n_r^{\\rm{post}}}}; \\frac{\\boldsymbol{\mathbf{\mu_r}} \circ \\boldsymbol{\mathbf{\pi_r}}}{\\boldsymbol{\mathbf{\mu_r}} \cdot \\boldsymbol{\mathbf{\pi_r}}} + \\boldsymbol{\mathbf{\\rho_r}} - \\boldsymbol{\mathbf{\delta_r}}\\right)

        \Pr\left(\\boldsymbol{\mathbf{n_r^{\\rm{err,pre}}}} \mid \\boldsymbol{\mathbf{\epsilon_r}}\\right) 
        &= \mathrm{Multinomial}\left(\\boldsymbol{\mathbf{n_r^{\\rm{err,pre}}}}; \\boldsymbol{\mathbf{\epsilon_r}}\\right)

        \Pr\left(\\boldsymbol{\mathbf{n_r^{\\rm{err,post}}}} \mid \\boldsymbol{\mathbf{\\rho_r}}\\right) 
        &= \mathrm{Multinomial}\left(\\boldsymbol{\mathbf{n_r^{\\rm{err,post}}}}; \\boldsymbol{\mathbf{\\rho_r}}\\right)

    where :math:`\\boldsymbol{\mathbf{\delta_r}}` is a vector that is 
    zero except for a one at the element corresponding to the wildtype.

    The MCMC tries to guarantee convergence via the parameters specified
    by `r_max`, `neff_min`, `nchains`, `niter`, `increasefactor`, 
    and `increasetries`.

    Args:
        `characterlist` (list)
           List of valid characters (e.g., codons, amino acids, nts).
        `wtchar` (str)
            Wildtype character at the site.
        `error_model` (object)
            Specifies how errors are estimated. Can be:
              - *none* : no errors (:math:`\\boldsymbol{\mathbf{\epsilon_r}} = \\boldsymbol{\mathbf{\\rho_r}} = \mathbf{0}`)
              - *same* : same error rates  pre and post-selection 
                (:math:`\\boldsymbol{\mathbf{\epsilon_r}} = \\boldsymbol{\mathbf{\\rho_r}}`).
              - *different* : different error rates pre- and post-selection
                (:math:`\\boldsymbol{\mathbf{\epsilon_r}} \\ne \\boldsymbol{\mathbf{\\rho_r}}`).
            Alternatively, can be instance of `StanModelNoneErr`, 
            `StanModelSameErr`, or `StanModelDifferentErr`, which implement
            the three strings described above. If you pass mode instances,
            then they will not need to be compiled, which saves time if
            calling this function repeatedly.
        `counts` (dict)
            Gives deep sequencing counts. Each string key should specify
            a dict keyed by all characters in `characterlist` with
            the values giving integer counts for that character. The keys:
                - `nrpre`: :math:`\\boldsymbol{\mathbf{n_r^{\\rm{pre}}}}`
                - `nrpost`: :math:`\\boldsymbol{\mathbf{n_r^{\\rm{post}}}}`
                - *nrerr*: required if `error_model` is `same`, specifies
                  :math:`\\boldsymbol{\mathbf{n_r^{\\rm{err,pre}}}} = \\boldsymbol{\mathbf{n_r^{\\rm{err,post}}}}`
                - `nrerrpre` and `nrerrpost`: required if `error_model` is
                  `different`, specify 
                  :math:`\\boldsymbol{\mathbf{n_r^{\\rm{err,pre}}}}` and
                  :math:`\\boldsymbol{\mathbf{n_r^{\\rm{err,post}}}}`.

        `priors` (dict)
            Specifies parameter vectors for Dirichlet priors. Each string
            key should specify a dict keyed by all characters and values
            giving the prior for that character. Keys:
                - `pir_prior_params`: :math:`\\boldsymbol{\mathbf{a_{\pi,r}}}`
                - `mur_prior_params`: :math:`\\boldsymbol{\mathbf{a_{\mu,r}}}`
                - `epsilonr_prior_params`: only required if `error_model` is 
                  `same` or `different`, specifies 
                  :math:`\\boldsymbol{\mathbf{a_{\epsilon,r}}}`
                - `rhor_prior_params` : only required if `error_model` is 
                  `different`, specifies
                  :math:`\\boldsymbol{\mathbf{a_{\\rho,r}}}`
            Any values less than `PRIOR_MIN_VALUE` (a module-level constant)
            are automatically set to `PRIOR_MIN_VALUE`.
        `seed` (int)
            Random number seed for MCMC. 
        `n_jobs` (int)
            Number of CPUs to use, -1 means all available.
        `niter`, `increasetries`, `r_max`, `neff_min`, `nchains`, `increasefactor`
            Specify MCMC convergence. They all have reasonable defaults.
            The MCMC is considered to have converged if the mean Gelman-Rubin R
            statistic (http://www.jstor.org/stable/2246093) over all 
            :math:`\pi_{r,x}` values <= `r_max` and the mean effective sample
            size is >= `neff_min`. The MCMC first runs with `nchains` chains
            and `niter` iterations. If it fails to converge, it increases
            the iterations by a factor of `increasefactor` and tries again,
            and repeats until it converges or has tried `increasetries`
            times to increase the number of iterations. If the effective
            sample size exceeds 3 times `neff_min` then we allow
            R to be `1 + 1.5 (r_max - 1)`.

    Returns:
        The tuple `(converged, pi_means, pi_95credint, logstring)` where:
            - `converged` is `True` if the MCMC converged, `False` otherwise.
            - `pi_means` is dict keyed by characters in `characterlist`
               with the value giving the :math:`\pi_{r,a}`.
            - `pi_95credint` is dict keyed by characters, values are 2-tuples
               giving median-centered credible interval for :math:`\pi_{r,a}`.
            - `logstring` is a string describing MCMC run and convergence.
    """
    logstring = ['\tBeginning MCMC at %s' % time.asctime()]
    numpy.random.seed(seed)
    assert nchains >= 2, "nchains must be at least two"
    assert niter >= 100, "niter must be at least 100"
    assert len(characterlist) == len(set(characterlist)), "characterlist contains a duplicate character:\n%s" % str(characterlist)
    assert wtchar in characterlist, "wtchar %s is not in characterlist %s" % (wtchar, str(characterlist))
    data = {'Nchar':len(characterlist), 'iwtchar':characterlist.index(wtchar) + 1}
    data['nrpre'] = [counts['nrpre'][char] for char in characterlist]
    data['nrpost'] = [counts['nrpost'][char] for char in characterlist]
    data['pir_prior_params'] = [max(PRIOR_MIN_VALUE, priors['pir_prior_params'][char]) for char in characterlist]
    data['mur_prior_params'] = [max(PRIOR_MIN_VALUE, priors['mur_prior_params'][char]) for char in characterlist]
    if error_model == 'none':
        sm = StanModelNoneErr().model
    elif isinstance(error_model, StanModelNoneErr):
        sm = error_model.model
        error_model = 'none'
    elif error_model == 'same' or isinstance(error_model, StanModelSameErr):
        if error_model == 'same':
            sm = StanModelSameErr().model
        else:
            sm = error_model.model
            error_model = 'same'
        data['nrerr'] = [counts['nrerr'][char] for char in characterlist]
        data['epsilonr_prior_params'] = [max(PRIOR_MIN_VALUE, priors['epsilonr_prior_params'][char]) for char in characterlist]
    elif error_model == 'different' or isinstance(error_model, StanModelDifferentErr):
        if error_model == 'different':
            sm = StanModelDifferentErr().model
        else:
            sm = error_model.model
            error_model = 'different'
        data['nrerrpre'] = [counts['nrerrpre'][char] for char in characterlist]
        data['nrerrpost'] = [counts['nrerrpost'][char] for char in characterlist]
        data['epsilonr_prior_params'] = [max(PRIOR_MIN_VALUE, priors['epsilonr_prior_params'][char]) for char in characterlist]
        data['rhor_prior_params'] = [max(PRIOR_MIN_VALUE, priors['rhor_prior_params'][char]) for char in characterlist]
    else:
        raise ValueError("Invalid error_model of %s" % error_model)
    ntry = 0
    while True: # run until converged or tries exhausted
        init = _InitialValuePreferences(error_model, nchains, characterlist.index(wtchar), len(characterlist))
        # this next command unfortunately creates a lot of output that I have not been able to re-direct...
        fit = sm.sampling(data=data, iter=niter, chains=nchains, seed=seed, n_jobs=n_jobs, refresh=-1, init=init)
        # extract output
        fitsummary = fit.summary()
        rownames = list(fitsummary['summary_rownames'])
        colnames = list(fitsummary['summary_colnames'])
        summary = fitsummary['summary']
        char_row_indices = dict([(char, rownames.index('pir[%d]' % characterlist.index(char))) for char in characterlist]) 
        rindex = colnames.index('Rhat')
        neffindex = colnames.index('n_eff')
        rlist = [summary[char_row_indices[char]][rindex] for char in characterlist]
        rhat_is_nan = [rhat for rhat in rlist if math.isnan(rhat)]
        rlist = [rhat for rhat in rlist if not math.isnan(rhat)]
        nefflist = [summary[char_row_indices[char]][neffindex] for char in characterlist]
        neffmean = sum(nefflist) / float(len(nefflist))
        if not rlist:
            assert len(rhat_is_nan) == len(characterlist)
            rmean = None
            logstring.append('\tAfter %d MCMC chains each of %d steps, mean R = nan and mean Neff = %d' % (nchains, niter, neffmean))
        else:
            rmean = sum(rlist) / float(len(rlist))
            logstring.append('\tAfter %d MCMC chains each of %d steps, mean R = %.2f and mean Neff = %d' % (nchains, niter, rmean, neffmean))
        if rhat_is_nan:
            logstring.append('\t\tThere are %d characters where R is nan' % len(rhat_is_nan))
        # allow convergence with extra stringent criteria when many sites have Rhat being nan
        # this is because pystan appears to sometimes give Rhat of nan for sites with low preference
        # turns out the pystan Rhat values of nan are a bug according to pystan developers 
        if (len(rhat_is_nan) < 0.25 * len(characterlist) and rmean != None and rmean <= r_max and neffmean >= neff_min) or (neffmean >= 3.0 * neff_min and ((rmean == None) or (rmean != None and rmean <= 1.0 + 1.5 * (r_max - 1.0)))):
            # converged
            logstring.append('\tMCMC is deemed to have converged at %s.' % time.asctime())
            meanindex = colnames.index('mean')
            lower95index = colnames.index('2.5%')
            upper95index = colnames.index('97.5%')
            pi_means = dict([(char, summary[char_row_indices[char]][meanindex]) for char in characterlist])
            pi_95credint = dict([(char, (summary[char_row_indices[char]][lower95index], summary[char_row_indices[char]][upper95index])) for char in characterlist])
            return (True, pi_means, pi_95credint, '\n'.join(logstring))
        else:
            # failed to converge
            if ntry < increasetries:
                ntry += 1
                niter = int(niter * increasefactor)
                logstring.append("\tMCMC failed to converge. Doing retry %d with %d iterations per chain." % (ntry, niter))
            else:
                with open('_no_converge_prefs_debugging.pickle', 'w') as f_debug:
                    pickle.dump((counts, init, fitsummary), f_debug)
                logstring.append("\tMCMC FAILED to converge after all attempts at %s." % time.asctime())
                meanindex = colnames.index('mean')
                lower95index = colnames.index('2.5%')
                upper95index = colnames.index('97.5%')
                pi_means = dict([(char, summary[char_row_indices[char]][meanindex]) for char in characterlist])
                pi_95credint = dict([(char, (summary[char_row_indices[char]][lower95index], summary[char_row_indices[char]][upper95index])) for char in characterlist])
                return (False, pi_means, pi_95credint, '\n'.join(logstring))




if __name__ == '__main__':
    import doctest
    doctest.testmod()
