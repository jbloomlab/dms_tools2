"""
========================
prefs
========================

Performs operations related to estimating site-specific amino-acid
preferences.

Uses `pystan <https://pystan.readthedocs.io/en/latest>`_
to perform MCMC for Bayesian inferences.
"""


import time
import math
import tempfile
import pickle
import random
import natsort
import numpy
import numpy.random
import pandas
import Bio.SeqIO
import pystan
import dms_tools2

#: minimum value for Dirichlet prior elements
PRIOR_MIN_VALUE = 1.0e-7 


class StanModelNoneErr(object):
    """``pystan`` model when `error_model` is `none`.
    
    For use by inferSitePrefs`."""
    def __init__(self, verbose=False):
        """Compile ``pystan`` model.

        Args:
            `verbose` (bool)
                Set to `True` if you want verbose compilation.
        """
        self.pystancode =\
"""
data {{
    int<lower=1> Nchar; // 64 for codons, 20 for amino acids, 4 for nucleotides
    int<lower=0> nrpre[Nchar]; // counts pre-selection
    int<lower=0> nrpost[Nchar]; // counts post-selection
    vector<lower={0:g}>[Nchar] pir_prior_params; // Dirichlet prior params
    vector<lower={0:g}>[Nchar] mur_prior_params; // Dirichlet prior params
}}
parameters {{
    simplex[Nchar] pir;
    simplex[Nchar] mur;
}}
transformed parameters {{
    simplex[Nchar] fr;
    fr = pir .* mur / dot_product(pir, mur);
}}
model {{
    pir ~ dirichlet(pir_prior_params);
    mur ~ dirichlet(mur_prior_params);
    nrpre ~ multinomial(mur);
    nrpost ~ multinomial(fr);
}}
""".format(PRIOR_MIN_VALUE)
        self.model = pystan.StanModel(model_code=self.pystancode, 
                verbose=verbose)


class StanModelSameErr:
    """``pystan`` model when `error_model` is `same`.
    
    For use by inferSitePrefs`."""
    def __init__(self, verbose=False):
        """Compile ``pystan`` model.

        Args:
            `verbose` (bool)
                Set to `True` if you want verbose compilation.
        """
        self.pystancode =\
"""
data {{
    int<lower=1> Nchar; // 64 for codons, 20 for amino acids, 4 for nucleotides
    int<lower=1, upper=Nchar> iwtchar; // index of wildtype character in 1, ... numbering
    int<lower=0> nrpre[Nchar]; // counts pre-selection
    int<lower=0> nrpost[Nchar]; // counts post-selection
    int<lower=0> nrerr[Nchar]; // counts in error control
    vector<lower={0:g}>[Nchar] pir_prior_params; // Dirichlet prior params
    vector<lower={0:g}>[Nchar] mur_prior_params; // Dirichlet prior params
    vector<lower={0:g}>[Nchar] epsilonr_prior_params; // Dirichlet prior params
}}
transformed data {{
    simplex[Nchar] deltar;
    for (ichar in 1:Nchar) {{
        deltar[ichar] = 0.0;
    }}
    deltar[iwtchar] = 1.0;
}}
parameters {{
    simplex[Nchar] pir;
    simplex[Nchar] mur;
    simplex[Nchar] epsilonr;
}}
transformed parameters {{
    simplex[Nchar] fr_plus_err;
    simplex[Nchar] mur_plus_err;
    fr_plus_err = pir .* mur / dot_product(pir, mur) + epsilonr - deltar;
    mur_plus_err = mur + epsilonr - deltar;
}}
model {{
    pir ~ dirichlet(pir_prior_params);
    mur ~ dirichlet(mur_prior_params);
    epsilonr ~ dirichlet(epsilonr_prior_params);
    nrerr ~ multinomial(epsilonr);
    nrpre ~ multinomial(mur_plus_err);
    nrpost ~ multinomial(fr_plus_err);
}}
""".format(PRIOR_MIN_VALUE)
        self.model = pystan.StanModel(model_code=self.pystancode,
                verbose=verbose)


class StanModelDifferentErr:
    """``pystan`` model when `error_model` is `different`.
    
    For use by inferSitePrefs`."""
    def __init__(self, verbose=False):
        """Compile ``pystan`` model.

        Args:
            `verbose` (bool)
                Set to `True` if you want verbose compilation.
        """
        self.pystancode =\
"""
data {{
    int<lower=1> Nchar; // 64 for codons, 20 for amino acids, 4 for nucleotides
    int<lower=1, upper=Nchar> iwtchar; // index of wildtype character in 1, ... numbering
    int<lower=0> nrpre[Nchar]; // counts pre-selection
    int<lower=0> nrpost[Nchar]; // counts post-selection
    int<lower=0> nrerrpre[Nchar]; // counts in pre-selection error control
    int<lower=0> nrerrpost[Nchar]; // counts in post-selection error control
    vector<lower={0:g}>[Nchar] pir_prior_params; // Dirichlet prior params
    vector<lower={0:g}>[Nchar] mur_prior_params; // Dirichlet prior params
    vector<lower={0:g}>[Nchar] epsilonr_prior_params; // Dirichlet prior params
    vector<lower={0:g}>[Nchar] rhor_prior_params; // Dirichlet prior params
}}
transformed data {{
    simplex[Nchar] deltar;
    for (ichar in 1:Nchar) {{
        deltar[ichar] = 0.0;
    }}
    deltar[iwtchar] = 1.0;
}}
parameters {{
    simplex[Nchar] pir;
    simplex[Nchar] mur;
    simplex[Nchar] epsilonr;
    simplex[Nchar] rhor;
}}
transformed parameters {{
    simplex[Nchar] fr_plus_err;
    simplex[Nchar] mur_plus_err;
    fr_plus_err = pir .* mur / dot_product(pir, mur) + rhor - deltar;
    mur_plus_err = mur + epsilonr - deltar;
}}
model {{
    pir ~ dirichlet(pir_prior_params);
    mur ~ dirichlet(mur_prior_params);
    epsilonr ~ dirichlet(epsilonr_prior_params);
    rhor ~ dirichlet(rhor_prior_params);
    nrerrpre ~ multinomial(epsilonr);
    nrerrpost ~ multinomial(rhor);
    nrpre ~ multinomial(mur_plus_err);
    nrpost ~ multinomial(fr_plus_err);
}}
""".format(PRIOR_MIN_VALUE)
        self.model = pystan.StanModel(model_code=self.pystancode,
                verbose=verbose)


def _initialValuePrefs(error_model, nchains, iwtchar, nchars):
    """Gets valid initial values for ``pystan`` preference inference.

    Values initialized by ``pystan`` frequently have invalid values.
    This function will generate random values that are valid
    for initialization and return a list that can be passed to the 
    `StanModel` as the `init` argument.
    """
    initattempts = 10 # might fail for pathological random values
    nrescales = 10 # rescale non-wildtype values down this many times
    rescalefactor = 5.0 # rescale non-wildtype values down by this much
    deltar = numpy.zeros(nchars)
    deltar[iwtchar] = 1.0
    init = []
    concparams = [1.0 for ichar in range(nchars)]
    concparams[iwtchar] = 10.0 # make this element bigger
    if error_model == 'none':
        for chain in range(nchains):
            while True:
                chain_init = {
                        'pir':numpy.random.dirichlet(concparams),
                        'mur':numpy.random.dirichlet(concparams),
                        }
                if (all(chain_init['pir'] > PRIOR_MIN_VALUE) and 
                        all(chain_init['mur'] > PRIOR_MIN_VALUE)):
                    init.append(chain_init)
                    break
        return init
    elif error_model == 'same':
        posterrname = 'epsilonr'
    elif error_model == 'different':
        posterrname = 'rhor'
    else:
        raise ValueError("Invalid error_model {0}".format(error_model))
    for chain in range(nchains):
        for iattempt in range(initattempts):
            chain_init = {
                    'pir':numpy.random.dirichlet(concparams),
                    'mur':numpy.random.dirichlet(concparams),
                    'epsilonr':numpy.random.dirichlet(concparams),
                    }
            if error_model == 'different':
                chain_init['rhor'] = numpy.random.dirichlet(concparams)
            irescale = 0
            while irescale < nrescales and (any(
                    chain_init['pir'] * chain_init['mur'] / 
                    numpy.dot(chain_init['pir'], chain_init['mur']) 
                    + chain_init[posterrname] - deltar < PRIOR_MIN_VALUE)
                    or any(chain_init['mur'] + chain_init['epsilonr'] 
                    - deltar < PRIOR_MIN_VALUE)): 
                irescale += 1
                chain_init['epsilonr'] /= rescalefactor
                chain_init['epsilonr'][iwtchar] = 1.0 - sum(
                        [chain_init['epsilonr'][ichar] for ichar 
                        in range(nchars) if ichar != iwtchar])
                if error_model == 'different':
                    chain_init['rhor'] /= rescalefactor
                    chain_init['rhor'][iwtchar] = 1.0 - sum([
                            chain_init['rhor'][ichar] for ichar in 
                            range(nchars) if ichar != iwtchar])
            if (all(chain_init['pir'] * chain_init['mur'] / 
                    numpy.dot(chain_init['pir'], chain_init['mur']) + 
                    chain_init[posterrname] - deltar > PRIOR_MIN_VALUE) 
                    and all(chain_init['mur'] + chain_init['epsilonr'] - 
                    deltar > PRIOR_MIN_VALUE)):
                break
        else:
            raise ValueError("Failed to initialize")
        init.append(chain_init)
    return init


def inferPrefsByRatio(charlist, sites, wts, pre, post, errpre,
        errpost, pseudocount):
    """Site-specific preferences from normalized enrichment ratios.

    Calculates site-specific preference :math:`\pi_{r,a}` of each
    site :math:`r` for each character :math:`a` using re-normalized
    enrichment ratios. This function accomplishes the same goal as
    `inferSitePrefs`, but is much simpler and faster, although less
    statistically rigorous in principle.

    Specifically, this function calculates the preferences as

    .. math::

        \pi_{r,a} = \\frac{\phi_{r,a}}{\sum_{a'} \phi_{r,a'}}

    where :math:`\phi_{r,a}` is the enrichment ratio of character
    :math:`a` relative to wildtype after selection. 
    
    The actual definition of these enrichment ratios is somewhat
    complex due to the need to correct for errors and add a
    pseudocount. We have 4 sets of counts:

        - `pre`: pre-selection counts

        - `post`: post-selection counts

        - `errpre`: error-control for pre-selection counts

        - `errpost`: error-control for post-selection counts

    For each set of counts :math:`s`, let :math:`N_r^s` (e.g., 
    :math:`N_r^{pre}`) denote the total counts for this sample at
    site :math:`r`, and let :math:`n_{r,a}^{s}` (e.g., 
    :math:`n_{r,a}^{pre}`) denote the counts at that site for
    character :math:`a`.

    Because some of the counts are low, we add a pseudocount 
    :math:`P` to each observation. Importantly, we need to scale
    this pseudocount by the total depth for that sample at that site
    in order to avoid systematically estimating different frequencies
    purely as a function of sequencing depth, which will bias
    the preference estimates.
    Therefore, given the pseudocount value :math:`P` defined via the
    ``--pseudocount`` argument to :ref:`dms2_prefs`, we calculate the
    scaled pseudocount for sample :math:`s` and site :math:`r` as

    .. math::

       P_r^s = P \\times \\frac{N_r^s}{\min\limits_{s'} N_r^{s'}}.

    This equation makes the pseudocount :math:`P` for the lowest-depth
    sample, and scales up the pseodocounts for all other samples
    proportional to their depth.

    With this pseudocount, the
    estimated frequency of character :math:`a` at site :math:`r`
    in sample :math:`s` is 
    
    .. math::
    
        f_{r,a}^s = \\frac{n_{r,a}^s + P_r^s}{N_{r}^s + A \\times P_r^s}

    where :math:`A` is the number of characters in the alphabet (e.g.,
    20 for amino acids without stop codons).

    The **error-corrected** estimates of the frequency of character
    :math:`a` before and after selection are then

    .. math::

        f_{r,a}^{before} &= \max\left(\\frac{P_r^{pre}}{N_{r}^{pre} + A \\times P_r^{pre}},
        \; f_{r,a}^{pre} + \delta_{a,\\rm{wt}\left(r\\right)}
        - f_{r,a}^{errpre}\\right)

        f_{r,a}^{after} &= \max\left(\\frac{P_r^{post}}{N_{r}^{post} + A \\times P_r^{post}},
        \; f_{r,a}^{post} + \delta_{a,\\rm{wt}\left(r\\right)}
        - f_{r,a}^{errpost}\\right)

    where :math:`\delta_{a,\\rm{wt}\left(r\\right)}` is the 
    Kronecker-delta, equal to 1 if :math:`a` is the same as the
    wildtype character :math:`\\rm{wt}\left(r\\right)` at site
    :math:`r`, and 0 otherwise. We use the :math:`\max` operator
    to ensure that even if the error-control frequency exceeds
    the actual estimate, our estimated frequency never drops
    below the pseudocount over the depth of the uncorrected sample.

    Given these definitions, we then simply calculate the enrichment 
    ratios as 

    .. math::

        \phi_{r,a} = \\frac{\left(f_{r,a}^{after}\\right) / 
        \left(f_{r,\\rm{wt}\left(r\\right)}^{after}\\right)}
        {\left(f_{r,a}^{before}\\right) / 
        \left(f_{r,\\rm{wt}\left(r\\right)}^{before}\\right)}

    In the case where we are **not** using any error-controls, then
    we simply set 
    :math:`f_{r,\\rm{wt}}^{errpre} = f_{r,\\rm{wt}}^{errpost} 
    = \delta_{a,\\rm{wt}\left(r\\right)}`.

    Args:
        `charlist` (list)
            Valid characters (e.g., codons, amino acids, nts).
        `sites` (list)
            Sites to analyze.
        `wts` (list)
            `wts[r]` is the wildtype character at site `sites[r]`.
        `pre` (pandas.DataFrame)
            Gives pre-selection counts. Should have columns
            with names of 'site' and all characters in `charlist`.
            The rows give the counts of each character at that site.
        `post` (pandas.DataFrame)
            Like `pre` but for post-selection counts.
        `errpre` (`None` or pandas.DataFrame)
            Like `pre` but for pre-selection error-control counts,
            or `None` if there is no such control.
        `errpost` (`None` or pandas.DataFrame)
            Like `pre` but for post-selection error-control counts,
            or `None` if there is no such control.
        `pseudocount` (float or int)
            The pseudocount to add to each observation.

    Returns:
        A pandas.DataFrame holding the preferences. The columns of
        this dataframe are 'site' and all characters in `charlist`.
        For each site in `sites`, the rows give the preference
        for that character.
    """
    assert len(wts) == len(sites) > 0
    assert all([wt in charlist for wt in wts]), "invalid char in wts"
    assert pseudocount > 0, "pseudocount must be greater than zero"

    dfs = {'pre':pre.copy(),
           'post':post.copy()
           }
    if errpre is not None:
        dfs['errpre'] = errpre.copy()
    if errpost is not None:
        dfs['errpost'] = errpost.copy()

    # compute total depth for each sample and sort by sites
    depths = []
    for stype in dfs.keys():
        df = dfs[stype]
        assert set(list(charlist) + ['site']) <= set(df.columns)
        assert set(sites) <= set(df['site'])
        df = df.query('site in @sites')
        assert len(df.index) == len(wts) == len(sites)
        df['site'] = pandas.Categorical(df['site'], sites)
        df = df.sort_values('site').reset_index(drop=True)
        df['Nr'] = df[charlist].sum(axis=1).astype('float')
        depths.append(df['Nr'])
        dfs[stype] = df
    mindepth = pandas.concat(depths, axis=1).min(axis=1)

    # calculate scaled pseudocounts
    pseudocounts = dict([(stype,
            pseudocount * (dfs[stype]['Nr'] / mindepth).fillna(1.0))
            for stype in dfs.keys()])

    # calculate frequencies
    for stype in dfs.keys():
        df = dfs[stype]
        for c in charlist:
            df['f{0}'.format(c)] = (df[c] + pseudocounts[stype]) / (
                    df['Nr'] + len(charlist) * pseudocounts[stype]) 
            df['delta{0}'.format(c)] = [int(wt == c) for wt in wts]
        dfs[stype] = df

    # fill values for errpre / errpost if they were None
    for stype in ['pre', 'post']:
        err = 'err{0}'.format(stype)
        if err not in dfs:
            dfs[err] = pandas.DataFrame(dict([('site', sites)] +
                    [('f{0}'.format(c), dfs[stype]['delta{0}'.format(c)])
                    for c in charlist]))

    # compute phi values
    fr = {}
    for (key, f, ferr, p) in [
            ('before', dfs['pre'], dfs['errpre'], pseudocounts['pre']),
            ('after', dfs['post'], dfs['errpost'], pseudocounts['post']),
            ]:
        fr[key] = {'wt':numpy.zeros(len(sites))}
        for c in charlist:
            fr[key][c] = numpy.maximum(
                    p / (f['Nr'].values + len(charlist) * p),
                    f['f{0}'.format(c)].values 
                        + f['delta{0}'.format(c)].values 
                        - ferr['f{0}'.format(c)].values
                    )
            fr[key]['wt'] += fr[key][c] * f['delta{0}'.format(c)].values
    phi = pandas.DataFrame({'site':sites})
    for c in charlist:
        phi[c] = ((fr['after'][c] / fr['after']['wt']) /
                  (fr['before'][c] / fr['before']['wt']))
    phi['denom'] = phi[charlist].sum(axis=1)

    # normalize phi values to get pi values
    prefs = phi[charlist].div(phi['denom'], axis=0)
    prefs.insert(0, 'site', sites)

    return prefs
    

def inferSitePrefs(charlist, wtchar, error_model, counts, 
        priors, seed=1, niter=10000, increasetries=5, n_jobs=1, 
        r_max=1.1, neff_min=100, nchains=4, increasefac=2):
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
    by `r_max`, `neff_min`, `nchains`, `niter`, `increasefac`, 
    and `increasetries`.

    Args:
        `charlist` (list)
           List of valid characters (e.g., codons, amino acids, nts).
        `wtchar` (str)
            Wildtype character at the site.
        `error_model` (str or object)
            Specifies how errors are estimated. Passing error model
            objects is faster if calling this function repeatedly
            as they will not need to be compiled. Can be:
              - The str `none` or instance of `StanModelNoneErr`:
                no errors (:math:`\\boldsymbol{\mathbf{\epsilon_r}} = \\boldsymbol{\mathbf{\\rho_r}} = \mathbf{0}`)
              - The str `same` or instance of `StanModelSameErr`: 
                same error rates  pre and post-selection 
                (:math:`\\boldsymbol{\mathbf{\epsilon_r}} = \\boldsymbol{\mathbf{\\rho_r}}`).
              - The str `different` or instance of `StanModelDifferentErr`: 
                different error rates pre- and post-selection
                (:math:`\\boldsymbol{\mathbf{\epsilon_r}} \\ne \\boldsymbol{\mathbf{\\rho_r}}`).
        `counts` (dict)
            Deep sequencing counts. Each string key should specify
            a dict keyed by all characters in `charlist` with
            the values giving integer counts for that character. The keys:
                - `pre`: :math:`\\boldsymbol{\mathbf{n_r^{\\rm{pre}}}}`
                - `post`: :math:`\\boldsymbol{\mathbf{n_r^{\\rm{post}}}}`
                - *err*: required if `error_model` is `same`, specifies
                  :math:`\\boldsymbol{\mathbf{n_r^{\\rm{err,pre}}}} = \\boldsymbol{\mathbf{n_r^{\\rm{err,post}}}}`
                - `errpre` and `errpost`: required if `error_model` is
                  `different`, specify 
                  :math:`\\boldsymbol{\mathbf{n_r^{\\rm{err,pre}}}}` and
                  :math:`\\boldsymbol{\mathbf{n_r^{\\rm{err,post}}}}`.
        `priors` (dict)
            Specifies parameter vectors for Dirichlet priors. Each string
            key should specify a dict keyed by all characters and values
            giving the prior for that character. Values less than 
            `PRIOR_MIN_VALUE` are set to `PRIOR_MIN_VALUE`. Keys are
                - `pir_prior_params`: :math:`\\boldsymbol{\mathbf{a_{\pi,r}}}`
                - `mur_prior_params`: :math:`\\boldsymbol{\mathbf{a_{\mu,r}}}`
                - `epsilonr_prior_params`: only required if `error_model` is 
                  `same` or `different`, specifies 
                  :math:`\\boldsymbol{\mathbf{a_{\epsilon,r}}}`
                - `rhor_prior_params` : only required if `error_model` is 
                  `different`, specifies
                  :math:`\\boldsymbol{\mathbf{a_{\\rho,r}}}`
        `seed` (int)
            Random number seed for MCMC. 
        `n_jobs` (int)
            Number of CPUs to use, -1 means all available.
        `niter`, `increasetries`, `r_max`, `neff_min`, `nchains`, `increasefac`
            Specify MCMC convergence. They all have reasonable defaults.
            The MCMC is considered to have converged if the mean Gelman-Rubin R
            statistic (http://www.jstor.org/stable/2246093) over all 
            :math:`\pi_{r,x}` values <= `r_max` and the mean effective sample
            size is >= `neff_min`. The MCMC first runs with `nchains` chains
            and `niter` iterations. If it fails to converge, it increases
            the iterations by a factor of `increasefac` and tries again,
            and repeats until it converges or has tried `increasetries`
            times to increase the number of iterations. If the effective
            sample size exceeds 3 times `neff_min` then we allow
            R to be `1 + 1.5 (r_max - 1)`.

    Returns:
        The tuple `(converged, pi_means, pi_95credint, logstring)` where:
            - `converged` is `True` if the MCMC converged, `False` otherwise.
            - `pi_means` is dict keyed by characters in `charlist`
               with the value giving the :math:`\pi_{r,a}`.
            - `pi_95credint` is dict keyed by characters, values are 2-tuples
               giving median-centered credible interval for :math:`\pi_{r,a}`.
            - `logstring` is a string describing MCMC run and convergence.
    """
    logstring = ['\tBeginning MCMC at %s' % time.asctime()]
    random.seed(seed)
    numpy.random.seed(seed)
    assert nchains >= 2, "nchains must be at least two"
    assert niter >= 100, "niter must be at least 100"
    assert len(charlist) == len(set(charlist))
    assert wtchar in charlist
    data = {'Nchar':len(charlist), 
            'iwtchar':charlist.index(wtchar) + 1,
            'nrpre':[counts['pre'][c] for c in charlist],
            'nrpost':[counts['post'][c] for c in charlist],
            'pir_prior_params':[max(PRIOR_MIN_VALUE, 
                    priors['pir_prior_params'][c]) for c in charlist],
            'mur_prior_params':[max(PRIOR_MIN_VALUE, 
                    priors['mur_prior_params'][c]) for c in charlist],
           }
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
        data['nrerr'] = [counts['err'][c] for c in charlist]
        data['epsilonr_prior_params'] = [max(PRIOR_MIN_VALUE, 
                priors['epsilonr_prior_params'][c]) for c in charlist]
    elif (error_model == 'different' or 
            isinstance(error_model, StanModelDifferentErr)):
        if error_model == 'different':
            sm = StanModelDifferentErr().model
        else:
            sm = error_model.model
            error_model = 'different'
        data['nrerrpre'] = [counts['errpre'][c] for c in charlist]
        data['nrerrpost'] = [counts['errpost'][c] for c in charlist]
        data['epsilonr_prior_params'] = [max(PRIOR_MIN_VALUE, 
                priors['epsilonr_prior_params'][c]) for c in charlist]
        data['rhor_prior_params'] = [max(PRIOR_MIN_VALUE, 
                priors['rhor_prior_params'][c]) for c in charlist]
    else:
        raise ValueError("Invalid error_model {0}".format(error_model))

    ntry = 0
    while True: # run until converged or tries exhausted
        init = _initialValuePrefs(error_model, nchains, 
                charlist.index(wtchar), len(charlist))
        fit = sm.sampling(data=data, iter=niter, chains=nchains, 
                seed=seed, n_jobs=n_jobs, refresh=-1, init=init)
        # extract output
        fitsummary = fit.summary()
        rownames = list(fitsummary['summary_rownames'])
        colnames = list(fitsummary['summary_colnames'])
        summary = fitsummary['summary']
        char_row_indices = dict([(c, rownames.index('pir[{0}]'.format(
                charlist.index(c)))) for c in charlist]) 
        rindex = colnames.index('Rhat')
        neffindex = colnames.index('n_eff')
        rlist = [summary[char_row_indices[c]][rindex] for c in charlist]
        rhat_is_nan = [rhat for rhat in rlist if math.isnan(rhat)]
        rlist = [rhat for rhat in rlist if not math.isnan(rhat)]
        nefflist = [summary[char_row_indices[c]][neffindex] for c in charlist]
        neffmean = sum(nefflist) / float(len(nefflist))
        if not rlist:
            assert len(rhat_is_nan) == len(charlist)
            rmean = None
            logstring.append('\tAfter {0} MCMC chains each of {1} steps, '
                    'mean R = nan and mean Neff = {2}'.format(
                    nchains, niter, neffmean))
        else:
            rmean = sum(rlist) / float(len(rlist))
            logstring.append('\tAfter {0} MCMC chains each of {1} steps, '
                    'mean R = {2} and mean Neff = {3}'.format(
                    nchains, niter, rmean, neffmean))
        if rhat_is_nan:
            logstring.append('\t\tThere are {0} characters where R is nan'
                    .format(len(rhat_is_nan)))
        # allow convergence with stringent criteria when Rhat is nan
        # pystan appears to give Rhat of nan for sites with low preference
        # pystan Rhat values of nan are a bug according to pystan developers
        if ((len(rhat_is_nan) < 0.25 * len(charlist) and 
                rmean != None and rmean <= r_max and neffmean >= neff_min) 
                or (neffmean >= 3.0 * neff_min and ((rmean == None) 
                or (rmean != None and rmean <= 1.0 + 1.5 * (r_max - 1.0))))):
            # converged
            logstring.append('\tMCMC converged at {0}.'.format(time.asctime()))
            meanindex = colnames.index('mean')
            lower95index = colnames.index('2.5%')
            upper95index = colnames.index('97.5%')
            pi_means = dict([(c, summary[char_row_indices[c]][meanindex]) 
                    for c in charlist])
            pi_95credint = dict([(c, 
                    (summary[char_row_indices[c]][lower95index], 
                    summary[char_row_indices[c]][upper95index])) 
                    for c in charlist])
            return (True, pi_means, pi_95credint, '\n'.join(logstring))
        else:
            # failed to converge
            if ntry < increasetries:
                ntry += 1
                niter = int(niter * increasefac)
                logstring.append("\tMCMC failed to converge. Doing retry "
                        "{0} with {1} iterations per chain.".format(
                        ntry, niter))
            else:
                with open('_no_converge_prefs_debug.pickle', 'wb') as f_debug:
                    pickle.dump((counts, init, fitsummary), f_debug)
                logstring.append("\tMCMC FAILED to converge after "
                        "all attempts at {0}.".format(time.asctime()))
                meanindex = colnames.index('mean')
                lower95index = colnames.index('2.5%')
                upper95index = colnames.index('97.5%')
                pi_means = dict([(c, summary[char_row_indices[c]][meanindex]) 
                        for c in charlist])
                pi_95credint = dict([(c, 
                        (summary[char_row_indices[c]][lower95index], 
                        summary[char_row_indices[c]][upper95index])) 
                        for c in charlist])
                return (False, pi_means, pi_95credint, '\n'.join(logstring))


def prefsToMutEffects(prefs, charlist):
    """Converts amino acid preferences to effects of specific mutations.

    If the preference of site :math:`r` for amino acid :math:`a` is
    :math:`\pi_{r,a}`, then the estimated effect (e.g., ratio of
    enrichment ratios) for mutating that site from :math:`x` to
    :math:`y` is :math:`\\frac{\pi_{r,y}}{\pi_{r,x}}`. Very small
    values indicate disfavored mutations, and very large values
    indicate highly favored mutations. The logarithm base 2 of the
    expected effects is also a useful measure -- negative values
    are disfavored mutations, positive values are favored ones.

    Args:
        `prefs` (pandas DataFrame)
            Preferences to analyze. The columns should be `site` and
            a column for every character in `charlist`.
        `charlist` (list)
            The list of characters that we are analyzing. For instance,
            `dms_tools2.AAS` for amino acids.

    Returns:
        A pandas Data Frame where the columns are:

            * `site`: the site

            * `initial`: the initial character (e.g., amino acid)

            * `final`: the final character after the mutation

            * `mutation`: mutation in the string form `A1G`

            * `effect`: the effect of the mutation

            * `log2effect`: the log2 of the effect of the mutation.

    >>> charlist = ['A', 'C', 'G']
    >>> prefs = pandas.DataFrame({
    ...         'site':[1, 2],
    ...         'A':[0.25, 0.25],
    ...         'C':[0.25, 0.5],
    ...         'G':[0.5, 0.25],
    ...         })
    >>> effects = prefsToMutEffects(prefs, charlist)
    >>> set(effects.columns) == {'site', 'initial', 'final',
    ...         'mutation', 'effect', 'log2effect'}
    True
    >>> numpy.allclose(effects[effects['initial'] == 'A']['effect'],
    ...         [1, 1, 2, 1, 2, 1])
    True
    >>> numpy.allclose(effects[effects['initial'] == 'C']['effect'],
    ...         [1, 1, 2, 0.5, 1, 0.5])
    True
    >>> numpy.allclose(effects['effect'], 2**effects['log2effect'])
    True
    """
    assert set(prefs.columns) <= set(['site'] + charlist)
    initial = prefs.melt(id_vars='site', value_vars=charlist,
                var_name='initial', value_name='initial_pref')
    final = initial.rename(columns=
            {'initial':'final', 'initial_pref':'final_pref'})
    effects = final.merge(initial, on='site')
    effects['effect'] = effects['final_pref'] / effects['initial_pref']
    effects['log2effect'] = numpy.log2(effects['effect'])
    effects['mutation'] = (effects['initial'] + 
            effects['site'].map(str) + effects['final'])
    return (effects
            .drop(['final_pref', 'initial_pref'], axis=1)
            .sort_values(['site', 'initial', 'final'])
            [['site', 'initial', 'final', 'mutation', 'effect', 'log2effect']]
            .reset_index(drop=True)
            )


def rescalePrefs(prefs, stringency):
    """Re-scale amino acid preferences by stringency parameter.

    If the initial preference of site :math:`r` for amino-acid 
    :math:`a` is :math:`\pi_{r,a}`, then the re-scaled preference is

    .. math::
    
        \\frac{\left(\pi_{r,a}\\right)^{\\beta}}{\sum_{a'} \left(\pi_{r,a'}\\right)}

    where :math:`\\beta` is the stringency parameter.

    Args:
        `prefs` (pandas.DataFrame)
            Columns are 'site' and then every character (e.g.,
            amino acid).
        `stringency` (float >= 0)
            The stringency parameter.

    Returns:
        A data frame in the same format as `prefs` but with the
        re-scaled preferences.

    >>> prefs = pandas.DataFrame(dict(
    ...         [('site', [1, 2]), ('A', [0.24, 0.03]), ('C', [0.04, 0.43])]
    ...         + [(aa, [0.04, 0.03]) for aa in dms_tools2.AAS if
    ...         aa not in ['A', 'C']]))
    >>> numpy.allclose(1, prefs.drop('site', axis=1).sum(axis=1))
    True
    >>> rescaled = rescalePrefs(prefs, 1.0)
    >>> all([numpy.allclose(prefs[c], rescaled[c]) for c in prefs.columns])
    True
    >>> rescaled2 = rescalePrefs(prefs, 2.0)
    >>> all(rescaled2['site'] == prefs['site'])
    True
    >>> numpy.allclose(rescaled2['A'], [0.6545, 0.0045], atol=1e-3)
    True
    >>> numpy.allclose(rescaled2['C'], [0.0182, 0.9153], atol=1e-3)
    True
    """
    assert set(prefs.drop('site', axis=1).select_dtypes(
            include=[numpy.number]).columns) == set(
            prefs.drop('site', axis=1).columns), ('Non-numeric '
            'columns other than "site"')

    chars = prefs.drop('site', axis=1).columns
    rescaled = prefs[chars].pow(stringency)
    rescaled = rescaled.divide(rescaled.sum(axis=1), axis='index')
    rescaled['site'] = prefs['site']
    return rescaled[prefs.columns]


def avgPrefs(prefsfiles):
    """Gets average of site-specific preferences.

    Args:
        `prefsfiles` (list)
            List of CSV files containing preferences, must all be
            for same sites and characters.

    Returns:
        A `pandas.DataFrame` containing the average of the
        preferences in `prefsfiles`. In this returned
        data frame, `site` is the index

    >>> tf1 = tempfile.NamedTemporaryFile
    >>> tf2 = tempfile.NamedTemporaryFile
    >>> with tf1(mode='w') as file1, tf2(mode='w') as file2:
    ...     x = file1.write('site,A,C,G,T\\n'
    ...                 '10,0.2,0.2,0.5,0.1\\n'
    ...                 '2a,0.3,0.3,0.3,0.1')
    ...     file1.flush()
    ...     x = file2.write('site,A,C,G,T\\n'
    ...                 '10,0.4,0.1,0.1,0.4\\n'
    ...                 '2a,0.3,0.4,0.1,0.2')
    ...     file2.flush()
    ...     avg = avgPrefs([file1.name, file2.name])
    >>> (avg['site'] == ['2a', '10']).all()
    True
    >>> numpy.allclose(avg['A'], [0.3, 0.3])
    True
    >>> numpy.allclose(avg['C'], [0.35, 0.15])
    True
    >>> numpy.allclose(avg['G'], [0.2, 0.3])
    True
    >>> numpy.allclose(avg['T'], [0.15, 0.25])
    True
    """
    assert len(prefsfiles) >= 1
    prefs = [pandas.read_csv(f, index_col='site').sort_index()
            for f in prefsfiles]

    # make sure all have the same columns in the same order
    cols = prefs[0].columns
    for i in range(len(prefs)):
        assert set(cols) == set(prefs[i].columns)
        prefs[i] = prefs[i][cols]

    avgprefs = pandas.concat(prefs).groupby('site').mean().reset_index()

    # natural sort by site: https://stackoverflow.com/a/29582718
    avgprefs = avgprefs.reindex(index=natsort.order_by_index(avgprefs.index,
            natsort.index_natsorted(avgprefs.site, signed=True)))

    return avgprefs


def prefsEntropy(prefs, charlist):
    """Calculate site entropy and number of effective characters.

    The site entropy :math:`h_r` for site :math:`r` is defined as
    :math:`h_r = -\sum_{x} \pi_{r,x} \log\left(\pi_{r,x}\\right)`
    where :math:`\pi_{r,x}` is the preference for character (e.g.,
    amino acid) :math:`x` at site :math:`r`, and the log is
    the natural logarithm.

    The number of effective characters at site :math:`r`
    is :math:`N_{\\rm{eff}, r} = \exp\left(h_r\\right)`.

    Args:
        `prefs` (pandas DataFrame)
            Data frame with the preferences. Should have a column
            named `site` and a column for each character in
            `charlist`. 
        `charlist` (list)
            List of the characters of interest, for example
            the amino acids.

    Returns:
        A copy of `prefs` that has additional columns named
        `entropy` and `neffective`, giving the site entropy and
        number of effective characters for each site. For
        the entropies, log is taken to base :math:`e`.

    >>> charlist = ['A', 'C', 'G', 'T']
    >>> prefs = pandas.DataFrame({
    ...         'site':[1, 2, 3],
    ...         'A':[0.25, 0.6, 1 / 3.],
    ...         'C':[0.25, 0.3, 1 / 3.],
    ...         'G':[0.25, 0.1, 1 / 3.],
    ...         'T':[0.25, 0.0, 0.0],
    ...         })
    >>> prefs_entropy = prefsEntropy(prefs, charlist)
    >>> (set(['entropy', 'neffective', 'site'] + charlist) == 
    ...         set(prefs_entropy.columns))
    True
    >>> h2 = -0.6 * math.log(0.6) - 0.3 * math.log(0.3) - 0.1 * math.log(0.1)
    >>> numpy.allclose(prefs_entropy['entropy'], [math.log(4), h2, math.log(3)])
    True
    >>> numpy.allclose(prefs_entropy['neffective'], [4, math.exp(h2), 3])
    True
    """
    assert 'site' in prefs.columns, "prefs does not have `site` column"
    assert set(charlist) < set(prefs.columns), "cols do not contain charlist"
    assert charlist, "no characters specified"
    assert numpy.allclose(prefs[charlist].sum(axis=1), 1, atol=1e-4),\
            "prefs do not sum to one"
    prefs_entropy = prefs.copy()
    prefs_entropy['entropy'] = 0
    for char in charlist:
        p_char = prefs_entropy[char].where(prefs_entropy[char] > 0, 1)
        prefs_entropy['entropy'] -= numpy.log(p_char) * p_char
    prefs_entropy['neffective'] = numpy.exp(prefs_entropy['entropy'])
    return prefs_entropy


def aafreqsFromAlignment(alignmentfile, codon_to_aa,
        ignore_gaps=True, ignore_stop=True):
    """Get amino-acid frequencies at each site in alignment.

    Args:
        `alignmentfile` (str)
            FASTA file with alignment of proteins or coding sequences.
        `codon_to_aa` (bool)
            If `True`, translate codon alignment to amino acids.
        `ignore_gaps` (bool)
            Ignore gaps when calculating frequencies.
        `ignore_stop` (bool)
            Ignore stop codons when calculating frequencies.

    Returns:
        A `pandas.DataFrame` with columns being `site` (1, 2, ...
        numbering) and other columns being amino acids and values
        giving frequencies in alignment.

    >>> with tempfile.NamedTemporaryFile(mode='w') as f:
    ...     x = f.write('>seq1\\n'
    ...                 'ATGGGGCAG\\n'
    ...                 '>seq2\\n'
    ...                 '---AGGCAG\\n'
    ...                 '>seq3\\n'
    ...                 'ATGTGACAG')
    ...     f.flush()
    ...     aafreqs = aafreqsFromAlignment(f.name, codon_to_aa=True)
    >>> aas_counts = ['M', 'G', 'R', 'Q']
    >>> aas_nocounts = [a for a in dms_tools2.AAS if a not in aas_counts]
    >>> (0 == aafreqs[aas_nocounts].values).all()
    True
    >>> expected_counts = pandas.DataFrame.from_items([
    ...         ('site', [1, 2, 3]), ('M', [1.0, 0.0, 0.0]),
    ...         ('G', [0.0, 0.5, 0]), ('R', [0.0, 0.5, 0.0]),
    ...         ('Q', [0.0, 0.0, 1.0])])
    >>> expected_counts.equals(aafreqs[['site'] + aas_counts])
    True
    """
    # read sequences
    seqs = [s.seq for s in Bio.SeqIO.parse(alignmentfile, 'fasta')]
    if codon_to_aa:
        seqs = [s.translate(gap='-', stop_symbol='*') for s in seqs]
    assert seqs, "No sequences"
    seqlen = len(seqs[0])
    assert seqlen, "sequences have no length"
    assert all([seqlen == len(s) for s in seqs]), "seqs not same length"

    # get character sets
    aas = dms_tools2.AAS.copy()
    skipchars = []
    if ignore_gaps:
        skipchars.append('-')
    else:
        aas.append('-')
    if ignore_stop:
        skipchars.append('*')
    else:
        aas.append('*')

    # tally amino-acid frequencies
    aafreqs = dict([(col, [0] * seqlen) for col in aas])
    aafreqs['site'] = list(range(1, seqlen + 1))
    for s in seqs:
        for (r, aa) in enumerate(s):
            if aa in skipchars:
                continue
            else:
                aafreqs[aa][r] += 1

    # convert to dataframe and change counts to freqs
    aafreqs = pandas.DataFrame(aafreqs)
    ncounts = aafreqs[aas].sum(axis=1).astype('float')
    for aa in aas:
        aafreqs[aa] = aafreqs[aa] / ncounts

    return aafreqs[['site'] + aas].fillna(0)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
