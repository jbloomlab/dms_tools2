"""
========================
compareprefs
========================

Performs operations related to comparing site-specific amino-acid
preferences among homologs.
"""

import dms_tools2
import math
import doctest

def ShannonJensenDivergence(p1, p2):
    """Returns the Shannon-Jensen divergence between *p1* and *p2*.

    This function returns the Shannon-Jensen divergence between the 
    probability distributions specified by *p1* and *p2*. The logarithms
    are taken to the base 2, meaning that the returned divergence
    will always be between 0 and 1. Requires ``numpy``, which is available
    if *PymcAvailable()* evaluates to *True*. Otherwise you will get an
    error.

    Args:
        `p1` and `p2` (*numpy.ndarray* objects)
            Two objects with elements that give the probability
            distributions for which we are computing the divergence.

    Returns:
        The Shannon-Jensen divergence (float)
    """
    assert 1 == len(p1.shape) == len(p2.shape), "p1 and p2 not both numpy.ndarrays of dimension 1"
    assert len(p1) == len(p2) > 0, "p1 and p2 not both arrays of same nonzero length"
    assert numpy.allclose(sum(p1), 1, atol=0.005), "p1 does not have entries summing to one: sum is %g, p1 is %s" % (sum(p1), p1)
    assert numpy.allclose(sum(p2), 1, atol=0.005), "p2 does not have entries summing to one: sum is %g, p2 is %s" % (sum(p2), p2)
    assert numpy.all(p1 >= 0), "p1 does not have all entries >= 0: p1 is %s" % p1
    assert numpy.all(p2 >= 0), "p2 does not have all entries >= 0: p2 is %s" % p2
    m = (p1 + p2) / 2.0
    d_p1_m = d_p2_m = 0.0
    for i in range(len(p1)):
        if p1[i]:
            d_p1_m += p1[i] * math.log(p1[i] / m[i], 2)
        if p2[i]:
            d_p2_m += p2[i] * math.log(p2[i] / m[i], 2)
    jsd = (d_p1_m + d_p2_m) / 2.0
    assert -1e10 <= jsd <= 1, "Shannon-Jensen divergence should be between zero and one, got value of %g" % jsd
    return jsd

def ComputeRMS(vector):
    """Computes the root mean square of all values in a vector
    
    Args:
        `vector` (list)
            A vector of values.
            
    Returns:
        The root mean square of the values in the input vector (float).

    >>> v = [1.2, 3.5, 6.8, 1.1]
    >>> round(ComputeRMS(v), 5)
    3.9096
    """
    RMS = math.sqrt( sum([math.pow(float(d),2) for d in vector]) / float(len(vector)) )
    return RMS

def MeasureDistanceBetweenVectors(v1, v2, distance_metric):
    """Computes the distance between two vectors of floats
    
    Args:
        `v1` and `v2` (lists)
            Vectors of floats to be compared
    
        `distance_metric` (string)
            Name of distance metric to be used. See below for options:
                * `half_sum_abs_diffs` : subtracts two vectors of floats, sums the
                absolute values of these differences, and then divides that value 
                by two so that the resulting value ranges between 0-1.
                * `normalized_RMSD` : computes the RMSD between two vectors of floats
                and multiplies this by the square root of 1/2 to normalize the
                resulting value such that it ranges between 0-1.
                * `Jensen_Shannon_distance` : the square root of the Jensen-Shannon
                divergence.
            
    Returns:
        The distance between the vectors (float)
        
    >>> v1 = [1.2, 3.5, 2.6]
    >>> v2 = [0.8, 4.0, 2.8]
    >>> round(MeasureDistanceBetweenVectors(v1, v2, 'half_sum_abs_diffs'),2)
    0.55
    >>> round(MeasureDistanceBetweenVectors(v1, v2, 'normalized_RMSD'), 2)
    0.27
    """
    assert len(v1) == len(v2)
    if distance_metric == 'half_sum_abs_diffs':
        distance = sum([abs(p1-p2) for (p1,p2) in zip(v1,v2)]) / 2.0
    elif distance_metric == 'normalized_RMSD':
        diffs = [p1-p2 for (p1,p2) in zip(v1,v2)]
        RMSD = ComputeRMS(diffs)
        distance = RMSD * math.sqrt(0.5)
    elif distance_metric == 'Jensen_Shannon_distance':
        distance = math.sqrt(ShannonJensenDivergence(numpy.asarray(v1), numpy.asarray(v2)))
    else:
        raise ValueError('Did not recognize distance metric: {0}'.format(distance_metric))
    return distance

def CompareHomologPrefs(inputfiles1, inputfiles2, outfile, sites=False, distance_metric='half_sum_abs_diffs'):
    """Measures site-specific differences in amino-acid preferences between two proteins
    
    This function compares the amino-acid preferences of one protein to the amino-acid
    preferences of a second protein at a site-by-site level. Specifically, at each site
    it measures how different that site's amino-acid preferences are between proteins,
    returning a large distance if the preferences are highly discordant, and a small
    distance if the preferences are nearly identical. Crucially, this metric also accounts
    for experimental noise by also considering differences among replicate measurements
    of the same protein.
    
    Args:
        `inputfiles1` (list or series)
            CSV files of replicate amino-acid preferences for homolog 1 in format returned by
            ``dms2_prefs`` preferences files.
        
        `inputfiles2` (list or series)
            Same as `inputfiles1`, but for homolog 2
            
        `sites` (list)
            Sites to compare. If no list is provided, all sites will be compared. Note that the
            input files must have all the same sites if no list is provided.

        `outfile` (string)
            Name of a CSV file with values of RMSDcorrected, RMSDwithin, RMSDbetween for each site

        `distance_metric` (string)
            Distance metric used in the comparison. See the function *MeasureDistanceBetweenVectors*
            for the different options. The default is *half_sum_abs_diffs*.
    
    Returns:
        A CSV file with values of RMSDcorrected, RMSDwithin, RMSDbetween for each site
    """
    # Assert that the number of replicates is greater than one for both samples
    nreps1 = len(inputfiles1)
    nreps2 = len(inputfiles2)
    assert nreps1 > 1 and nreps2 > 1, 'Need at least one input file for each homolog'
    
    # Read in amino-acid preferences and store as dataframes in a dictionary
    # keyed by homolog and then by replicate. If `sites` is included as an input,
    # only read in those sites. Otherwise consider all sites.
    prefs = {}
    sites_in_each_dataset = []
    for (h, inputfiles) in [('h1', inputfiles1), ('h2', inputfiles2)]:
        prefs[h] = {}
        for (i, f) in enumerate(inputfiles, 1):
            if sites:
                sites = list(map(str, sites)) # convert list to a list of strings if necessary
                p = pandas.read_csv(f, index_col='site')
                assert all([site in p.index.values for site in sites]
                          ), "Not all of the listed sites are in the input file {0}".format(f)
                prefs[h][i] = p.loc[sites]
            else:
                prefs[h][i] = pandas.read_csv(f, index_col='site')
            sites_in_each_dataset.append(list(prefs[h][i].index.values))    
    
    # Make sure each dataset has the exact same set of sites
    sites = sites_in_each_dataset[0]
    assert all([sites == i for i in sites_in_each_dataset]
              ), "Not all input files have the exact same set of sites"
    
    # Make list of comparisons within and between homologs. I will keep track of each
    # comparison as a tupple of tupples. E.g., ((h1,r1), (h1,r2)) encodes the comparison
    # of replicates 1 and 2 for the first homolog. This method will come in handy when
    # keying the above dictionary of amino-acid preferences.
    hr1 = [('h1', r) for r in prefs['h1']] # each sample is a tupple: (homolog, replicate)
    hr2 = [('h2', r) for r in prefs['h2']]
    within_comparisons = list(itertools.combinations(hr1, 2)) + list(itertools.combinations(hr2, 2))
    all_comparisons = list(itertools.combinations(hr1 + hr2, 2))
    between_comparisons = [comparison for comparison in all_comparisons if comparison not in within_comparisons]
    assert len(set.intersection(set(within_comparisons), set(between_comparisons))) == 0
    
    # Compute site-specific distances at each site
    distances = dict([(key, []) for key in ['site', 'RMSD_within', 'RMSD_between', 'RMSD_corrected']])
    for site in sites:
        within_distances = []
        between_distances = []
        for comparison in all_comparisons:
            ((h1, r1), (h2, r2)) = comparison
            prefs1 = prefs[h1][r1].loc[site]
            prefs2 = prefs[h2][r2].loc[site]
            distance = MeasureDistanceBetweenVectors(prefs1, prefs2, distance_metric)
            if comparison in within_comparisons:
                within_distances.append(distance)
            elif comparison in between_comparisons:
                between_distances.append(distance)
            else:
                raise ValueError('Comparison not in lists of within or between comparisons: {0}'.format(comparison))
        assert len(within_comparisons) == len(within_distances)
        assert len(between_comparisons) == len(between_distances)
        RMSD_within = ComputeRMS(within_distances)
        RMSD_between = ComputeRMS(between_distances)
        RMSD_corrected = RMSD_between - RMSD_within
        distances['site'].append(site)
        distances['RMSD_within'].append(RMSD_within)
        distances['RMSD_between'].append(RMSD_between)
        distances['RMSD_corrected'].append(RMSD_corrected)
        
    # Write results to an outfile
    distances_df = pandas.DataFrame.from_dict(distances)
    distances_df.set_index('site', inplace=True)
    distances_df[['RMSD_corrected', 'RMSD_between', 'RMSD_within']].sort_values(by='RMSD_corrected', ascending=False).to_csv(outfile)
    
    pass # Does not return anything aside from the output file


if __name__ == '__main__':
    import doctest
    doctest.testmod()
