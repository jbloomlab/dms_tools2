"""
========================
compareprefs
========================

Compare site-specific amino-acid preferences among homologs.
"""

import math
import tempfile
import io
import functools
import itertools

import natsort
import numpy
import pandas

import dms_tools2


def divJensenShannon(p1, p2):
    """Jensen-Shannon divergence between two distributions.

    The logarithms are taken to base 2, so the result will be
    between 0 and 1.

    Args:
        `p1` and `p2` (array-like)
            The two distributions for which we compute divergence.

    Returns:
        The Jensen-Shannon divergence as a float.

    >>> p1 = [0.5, 0.2, 0.2, 0.1]
    >>> p2 = [0.4, 0.1, 0.3, 0.2]
    >>> p3 = [0.0, 0.2, 0.2, 0.6]
    >>> numpy.allclose(divJensenShannon(p1, p1), 0, atol=1e-5)
    True
    >>> numpy.allclose(divJensenShannon(p1, p2), 0.035789, atol=1e-5)
    True
    >>> numpy.allclose(divJensenShannon(p1, p3), 0.392914, atol=1e-5)
    True
    """
    p1 = numpy.array(p1)
    p2 = numpy.array(p2)

    def _kldiv(a, b):
        with numpy.errstate(all='ignore'):
            kl = numpy.sum([v for v in a * numpy.log2(a / b)
                    if not numpy.isnan(v)])
        return kl

    m = 0.5 * (p1 + p2)

    return 0.5 * (_kldiv(p1, m) +_kldiv(p2, m))


def computeRMS(v):
    """The root mean square (RMS) of a list of values
    
    Args:
        `v` (array-like)
            Values for which we compute the RMS.
            
    Returns:
        The RMS of the values.

    >>> v = [1.2, 3.5, 6.8, 1.1]
    >>> numpy.allclose(computeRMS(v), 3.9096, atol=1e-4)
    True
    """
    v = numpy.array(v)
    assert len(v) > 0
    return math.sqrt((v**2).sum() / len(v))


def prefDistance(pi1, pi2, distmetric):
    """Computes the distance between two arrays of preferences.
    
    Args:
        `pi1` and `pi2` (array-like)
            Two arrays of preferences.
        `distmetric` (string)
            Distance metric to use. Can be:
                - `half_sum_abs_diff`: half sum absolute value of difference
                - `JensenShannon`: square root of Jensen-Shannon divergence
            
    Returns:
        The distance between `pi1` and `pi2`.
        
    >>> pi1 = [0.5, 0.2, 0.3]
    >>> pi2 = [0.2, 0.4, 0.4]
    >>> numpy.allclose(prefDistance(pi1, pi1, 'half_sum_abs_diff'), 0)
    True
    >>> numpy.allclose(prefDistance(pi1, pi1, 'JensenShannon'), 0)
    True
    >>> numpy.allclose(prefDistance(pi1, pi2, 'half_sum_abs_diff'), 0.3)
    True
    >>> numpy.allclose(prefDistance(pi1, pi2, 'JensenShannon'), 0.2785483)
    True
    """
    pi1 = numpy.array(pi1)
    pi2 = numpy.array(pi2)
    assert len(pi1) == len(pi2)
    assert numpy.allclose(pi1.sum(), 1, atol=0.005)
    assert numpy.allclose(pi2.sum(), 1, atol=0.005)
    assert numpy.all(pi1 >= 0)
    assert numpy.all(pi2 >= 0)

    if distmetric == 'half_sum_abs_diff':
        dist = (numpy.fabs(pi1 - pi2)).sum() / 2.0

    elif distmetric == 'JensenShannon':
        dist = math.sqrt(divJensenShannon(pi1, pi2))

    else:
        raise ValueError('Invalid `distmetric` {0}'.format(distmetric))

    return dist


def comparePrefs(prefs1, prefs2, sites=None, distmetric='half_sum_abs_diff',
        chars=dms_tools2.AAS):
    """Compute error-corrected distance between two sets of preferences.
  
    Designed for the situation in which you have made replicate
    measurements of the amino-acid preferences for two protein
    homologs, and want to estimate the difference in preferences
    at each site while correcting for experimental error as
    quantified by the replicate measurements.

    The *distance* between each pair of replicates at each
    site is computed using `prefDistance` with `distmetric`.
    We then compute the RMS distance between all pairs
    for the same homolog to get `RMSDwithin`, and all pairs
    of different homologs to get `RMSDbetween`. We calculate
    `RMSDcorrected` as `RMSDbetween - RMSDwithin`.

    We also compute the mean (across replicates) preference
    for homolog 1 minus the mean for homolog 2, scaled so
    that the total height in each direction equals `RMSDcorrected`.
    These values are an error-corrected estimate of the difference
    in preference for each amino acid between homologs.

    Args:
        `prefs1` (list)
            Files giving replicate measurements of preferences for
            homolog 1 in the CSV format returned by ``dms2_prefs``.
        `prefs2` (list)
            Files giving measurements for homolog 2.
        `sites` (list or `None`)
            If `None`, compare all sites shared between the two
            homolog preference sites. Otherwise should be a list
            of the sites to compare.
        `distmetric` (string)
            Distance metric to use. Can be any valid option for
            the argument of the same name to `prefDistance`.
        `chars` (list)
            List of characters for which we analyze the preferences.
            For instance, all 20 amino acids.
    
    Returns:
        A `pandas.DataFrame` giving the distances at each site,
        as well as the replicate mean difference between 
        preferences for homolog 1 minus homolog 2 for each amino
        acid at each site scaled to height of `RMSDcorrected`
        in each direction.

    Example calculation for two character sequences and two
    replicates for each homolog:

    >>> TF = functools.partial(tempfile.NamedTemporaryFile, mode='w')
    >>> with TF() as p1_1, TF() as p1_2, TF() as p2_1, TF() as p2_2:
    ...     n = p1_1.write('''site,    A,    C
    ...                          1,  0.8,  0.2
    ...                          2,  0.3,  0.7'''.replace(' ', ''))
    ...     p1_1.flush()
    ...     n = p1_2.write('''site,    A,    C
    ...                          1,  0.8,  0.2
    ...                          2,  0.4,  0.6'''.replace(' ', ''))
    ...     p1_2.flush()
    ...     n = p2_1.write('''site,    A,    C
    ...                          2,  0.4,  0.6
    ...                          1,  0.6,  0.4'''.replace(' ', ''))
    ...     p2_1.flush()
    ...     n = p2_2.write('''site,    A,    C
    ...                          1,  0.6,  0.4
    ...                         1a,  0.4,  0.6
    ...                          2,  0.5,  0.5'''.replace(' ', ''))
    ...     p2_2.flush()
    ...     diffs = comparePrefs([p1_1.name, p1_2.name],
    ...                          [p2_1.name, p2_2.name],
    ...                          chars=['A', 'C'])
    >>> print(diffs.to_string(float_format=lambda x: '{0:.2f}'.format(x)))
      site  RMSDcorrected  RMSDbetween  RMSDwithin     A     C
    0    1           0.20         0.20        0.00  0.20 -0.20
    1    2           0.02         0.12        0.10 -0.02  0.02
    """

    assert len(prefs1) > 1, "provide prefs for multiple replicates"
    assert len(prefs2) > 1, "provide prefs for multiple replicates"
  
    # read in all preferences
    prefs = []
    expectcols = ['site'] + chars
    for (homolog, homologprefs) in enumerate([prefs1, prefs2], 1):
        for (rep, repprefs) in enumerate(homologprefs, 1):
            iprefs = pandas.read_csv(repprefs)
            iprefs['site'] = iprefs['site'].astype('str')
            assert set(iprefs.columns) <= set(expectcols), \
                    "{0} missing expected columns".format(repprefs)
            prefs.append(iprefs[expectcols]
                               .assign(homolog=homolog, replicate=rep)
                               )

    # get only desired sites
    if sites is None:
        # use sites shared among all preference sets
        sites = list(set.intersection(*[set(p['site']) for p in prefs]))
    assert isinstance(sites, list) and len(sites), "no `sites` to analyze"
    sites = natsort.natsorted(list(map(str, sites)), signed=True)

    # merge preferences for desired sites
    assert all([set(p['site']) >= set(sites) for p in prefs]),\
            "not all prefs have all sites"
    prefs = [p[p['site'].isin(sites)] for p in prefs]
    prefs = pandas.concat(prefs)
    prefs['site'] = pandas.Categorical(prefs['site'], sites)
    prefs = prefs.sort_values('site').set_index('site')

    # compute RMSDs
    dists = {'within':[], 'between':[]}
    for ((hi, repi), (hj, repj)) in itertools.combinations(
            [(h, rep) for h in [1, 2] for rep in
            prefs.query('homolog == @h')['replicate'].unique()], 2):
        prefsi = (prefs.query('homolog == @hi and replicate == @repi')
                  [chars])
        prefsj = (prefs.query('homolog == @hj and replicate == @repj')
                  [chars])
        assert prefsi.index.equals(prefsj.index)
        disttype = {True:'within', False:'between'}[hi == hj]
        dists[disttype].append(
                [prefDistance(prefsi.loc[r], prefsj.loc[r], distmetric)
                for r in sites]
                )
    for (disttype, dist) in dists.items():
        distseries = (pandas.DataFrame(dist, columns=sites)
                      .transpose()
                      .apply(computeRMS, axis=1)
                      )
        prefs['RMSD' + disttype] = distseries
    prefs['RMSDcorrected'] = prefs['RMSDbetween'] - prefs['RMSDwithin']
    rmsds = ['RMSDcorrected', 'RMSDbetween', 'RMSDwithin']

    # compute RMSDcorrected-scaled diff between homologs for each pref
    prefmeans = {}
    for homolog in [1, 2]:
        prefmeans[homolog] = (prefs.reset_index()
                              .query('homolog == @homolog')
                              .groupby('site')
                              [chars]
                              .mean()
                              )
    prefs = prefs[~prefs.index.duplicated(keep='first')][rmsds]
    dprefs = prefmeans[1] - prefmeans[2]
    # normalize so sums to one in each direction
    dprefs = dprefs.div(dprefs.abs().sum(axis=1), axis=0).mul(2).fillna(0)
    dprefs = dprefs.mul(prefs['RMSDcorrected'], axis=0)
    prefs = prefs.join(dprefs)

    return prefs[rmsds + chars].reset_index()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
