"""
===========
fracsurvive
===========

Performs operations related to estimating the fraction of each
mutation that survives a selection.
"""

import os
import tempfile
import numpy
import pandas
from dms_tools2 import CODONS, CODON_TO_AA


def computeMutFracSurvive(libfracsurvive, sel, mock, countcharacters,
        pseudocount, translate_to_aa, err=None, mincount=0):
    """Compute fraction surviving for each mutation.

    Args:
        `libfracsurvive` (float)
            Overall fraction of selected library that survives
            relative to mock-selected. Should be >= 0 and <= 1.
        `sel` (pandas.DataFrame)
            Counts for selected sample. Columns should be
            `site`, `wildtype`, and every character in
            `countcharacters`.
        `mock` (pandas.DataFrame)
            Like `sel` but counts for mock-selected sample.
        `countcharacters` (list)
            List of all characters (e.g., codons).
        `pseudocount` (int or float > 0)
            Pseudocount to add to counts.
        `translate_to_aa` (bool)
            Should be `True` if counts are for codons and we are
            estimating mutfracsurvive for amino acids, `False` otherwise.
        `err` (pandas.DataFrame or `None`)
            Optional error-control counts, in same format as `sel`.
        `mincount` (int >= 0)
            Report as `NaN` the mutfracsurvive for any mutation where 
            neither `sel` nor `mock` has at least this many counts.

    Returns:
        A `pandas.DataFrame` with the fraction surviving for each mutation.
        Columns are `site`, `wildtype`, `mutation`, `mutfracsurvive`.

    >>> countchars = ['A', 'C', 'G', 'T']
    >>> libfracsurvive = 0.1
    >>> pseudocount = 5
    >>> mock = pandas.DataFrame.from_records(
    ...         [(1, 'A', 95, 95, 95, 95), (2, 'C', 195, 195, 95, 95)],
    ...         columns=['site', 'wildtype', 'A', 'C', 'G', 'T'])
    >>> sel = pandas.DataFrame.from_records(
    ...         [(1, 'A', 390, 90, 90, 190), (2, 'C', 390, 190, 390, 190)],
    ...         columns=['site', 'wildtype', 'A', 'C', 'G', 'T'])
    >>> mutfracsurvive = computeMutFracSurvive(libfracsurvive, sel, mock,
    ...         countchars, pseudocount, False)
    >>> {'site', 'wildtype', 'mutation', 'mutfracsurvive'} == set(
    ...         mutfracsurvive.columns)
    True
    >>> mutfracsurvive = mutfracsurvive.sort_values(['site', 'mutation'])
    >>> all(mutfracsurvive['site'] == [1, 1, 1, 1, 2, 2, 2, 2])
    True
    >>> all(mutfracsurvive['mutation'] == countchars + countchars)
    True
    >>> numpy.allclose(mutfracsurvive.query('site == 1')['mutfracsurvive'],
    ...         [0.2, 0.05, 0.05, 0.1])
    True
    >>> numpy.allclose(mutfracsurvive.query('site == 2')['mutfracsurvive'],
    ...         [0.1, 0.05, 0.2, 0.1])
    True
    """
    assert pseudocount > 0

    assert 0 <= libfracsurvive <= 1

    expectedcols = set(['site', 'wildtype'] + countcharacters)
    assert set(sel.columns) == expectedcols, "Invalid columns for sel"
    assert set(mock.columns) == expectedcols, "Invalid columns for sel"

    sel = sel.sort_values('site')
    mock = mock.sort_values('site')
    assert all(sel['site'] == mock['site']), "Inconsistent sites"
    assert all(sel['wildtype'] == mock['wildtype']), "Inconsistent wildtype"

    if err is not None:
        assert (set(err.columns) == expectedcols), "Invalid columns for err"
        err = err.sort_values('site')
        assert all(err['site'] == sel['site']), "Inconsistent sites"
        assert all(err['wildtype'] == sel['wildtype']), "Inconsistent sites"

    # make melted dataframe for calculations
    m = None    
    for (df, name) in [(sel, 'sel'), (mock, 'mock'), (err, 'err')]:
        if (name == 'err') and (err is None):
            continue
        m_name = (df.melt(id_vars=['site', 'wildtype'],
                          var_name='mutation', value_name='n')
                    .sort_values(['site', 'mutation'])
                    .reset_index(drop=True)
                    )
        m_name['N'] = m_name.groupby('site').transform('sum')['n']
        if m is None:
            m = m_name[['site', 'wildtype', 'mutation']].copy()
        assert all([all(m[c] == m_name[c]) for c in 
                ['site', 'wildtype', 'mutation']])
        m['n{0}'.format(name)] = m_name['n'].astype('float')
        m['N{0}'.format(name)] = m_name['N'].astype('float')

    # error correction 
    if err is not None:
        m['epsilon'] = m['nerr'] / m['Nerr']
        wtmask = m['mutation'] == m['wildtype']
        assert all(m[wtmask] > 0), "err counts of 0 for wildtype"
        for name in ['sel', 'mock']:
            ncol = 'n{0}'.format(name)
            Ncol = 'N{0}'.format(name)
            # follow here to set wildtype values without zero division
            wtval = pandas.Series(0, wtmask.index)
            wtval.loc[wtmask] = m.loc[wtmask, ncol] / m.loc[wtmask, 'epsilon']
            m[ncol] = numpy.maximum(0, 
                    (m[Ncol] * (m[ncol] / m[Ncol] - m['epsilon']))
                    ).where(~wtmask, wtval)
            m[Ncol] = m.groupby('site').transform('sum')[ncol]
        m = m.reset_index()
   
    # convert codon to amino acid counts
    if translate_to_aa:
        assert set(countcharacters) == set(CODONS),\
                "translate_to_aa specified, but not using codons"
        m = (m.replace({'wildtype':CODON_TO_AA, 'mutation':CODON_TO_AA})
              .groupby(['site', 'wildtype', 'mutation', 'Nsel', 'Nmock'])
              .sum()
              .reset_index()
              )

    # add pseudocounts
    m['nselP'] = (m['nsel'] + pseudocount * numpy.maximum(1,
            m['Nsel'] / m['Nmock']))
    m['nmockP'] = (m['nmock'] + pseudocount * numpy.maximum(1,
            m['Nmock'] / m['Nsel']))
    m['NselP'] = (m['Nsel'] + pseudocount * numpy.maximum(1,
            m['Nsel'] / m['Nmock']) * len(countcharacters))
    m['NmockP'] = (m['Nmock'] + pseudocount * numpy.maximum(1,
            m['Nmock'] / m['Nsel']) * len(countcharacters))

    # compute mutfracsurvive
    m['mutfracsurvive'] = (libfracsurvive *
            (m['nselP'] / m['NselP']) / (m['nmockP'] / m['NmockP'])
            ).where((m['nsel'] >= mincount) | (m['nmock'] >= mincount),
                    numpy.nan)

    return m[['site', 'wildtype', 'mutation', 'mutfracsurvive']]


def mutToSiteFracSurvive(mutfracsurvive):
    """Computes sitefracsurvive from mutfracsurvive.

    Args:
        `mutfracsurvive` (pandas.DataFrame)
            Dataframe with mutfracsurvive as from `computeMutFracSurvive`

    Returns:
        The dataframe `sitefracsurvive`, which has the following columns:
            - `site`
            - `sitefracsurvive`: sum of mutfracsurvive at site
            - `maxfracsurvive`: maximum mutfracsurvive at site

    >>> mutfracsurvive = (pandas.DataFrame({
    ...         'site':[1, 2, 3, 4],
    ...         'wildtype':['A', 'C', 'C', 'A'],
    ...         'A':[numpy.nan, 0.9, 0.1, numpy.nan],
    ...         'C':[0.2, numpy.nan, numpy.nan, 0.3],
    ...         'G':[0.8, 0.1, 0.2, 0.0],
    ...         'T':[0.2, 0.0, 0.1, 0.4],
    ...         })
    ...         .melt(id_vars=['site', 'wildtype'],
    ...               var_name='mutation', value_name='mutfracsurvive')
    ...         .reset_index(drop=True)
    ...         )
    >>> sitefracsurvive = mutToSiteFracSurvive(mutfracsurvive)
    >>> all(sitefracsurvive.columns == ['site', 'sitefracsurvive', 
    ...         'maxfracsurvive'])
    True
    >>> all(sitefracsurvive['site'] == [1, 2, 3, 4])
    True
    >>> numpy.allclose(sitefracsurvive['sitefracsurvive'],
    ...         [1.2, 1.0, 0.4, 0.7])
    True
    >>> numpy.allclose(sitefracsurvive['maxfracsurvive'],
    ...         [0.8, 0.9, 0.2, 0.4])
    True
    """
    sitefracsurvive = (
            mutfracsurvive
            .pivot(index='site', columns='mutation', values='mutfracsurvive')
            .fillna(0)
            .assign(sitefracsurvive=lambda x: x.abs().sum(axis=1),
                    maxfracsurvive=lambda x: x.max(axis=1),
                    )
            .reset_index()
            [['site', 'sitefracsurvive', 'maxfracsurvive']]
            )
    return sitefracsurvive


def avgMutFracSurvive(mutfracsurvivefiles, avgtype):
    """Gets mean or median mutation fraction surviving.

    Typically used to get an average across replicates.

    Args:
        `mutfracsurvivefiles` (list)
            List of CSV files with mutfracsurvivesel as returned by
            ``dms2_fracsurvive``.
        `avgtype` (str)
            Type of "average" to calculate. Possibilities:
                - `mean`
                - `median`

    Returns:
        A `pandas.DataFrame` containing the mean or median
        mutation fraction survive (`mutfracsurvive`).

    >>> tf = tempfile.NamedTemporaryFile
    >>> with tf(mode='w') as f1, tf(mode='w') as f2, tf(mode='w') as f3:
    ...     x = f1.write('site,wildtype,mutation,mutfracsurvive\\n'
    ...                  '156,G,S,0.9\\n'
    ...                  '157,K,D,0.1')
    ...     f1.flush()
    ...     x = f2.write('site,wildtype,mutation,mutfracsurvive\\n'
    ...                  '157,K,D,0.1\\n'
    ...                  '156,G,S,1.0')
    ...     f2.flush()
    ...     x = f3.write('site,wildtype,mutation,mutfracsurvive\\n'
    ...                  '157,K,D,0.4\\n'
    ...                  '156,G,S,0.5')
    ...     f3.flush()
    ...     mean = avgMutFracSurvive([f1.name, f2.name, f3.name],
    ...             'mean')
    ...     median = avgMutFracSurvive([f1.name, f2.name, f3.name],
    ...             'median')
    >>> (mean['site'] == [156, 157]).all()
    True
    >>> (median['site'] == [156, 157]).all()
    True
    >>> (mean['wildtype'] == ['G', 'K']).all()
    True
    >>> (median['wildtype'] == ['G', 'K']).all()
    True
    >>> (mean['mutation'] == ['S', 'D']).all()
    True
    >>> (median['mutation'] == ['S', 'D']).all()
    True
    >>> numpy.allclose(mean['mutfracsurvive'], [0.8, 0.2])
    True
    >>> numpy.allclose(median['mutfracsurvive'], [0.9, 0.1])
    True
    """
    fracsurvives = []
    cols = ['site', 'wildtype', 'mutation', 'mutfracsurvive']
    for f in mutfracsurvivefiles:
        fracsurvives.append(
                pandas.read_csv(f)
                [cols]
                .sort_values(['site', 'wildtype', 'mutation'])
                .reset_index()
                )
        assert all([(fracsurvives[0][c] == fracsurvives[-1][c]).all()
                for c in ['site', 'wildtype', 'mutation']]),\
                "files do not have same muts" 
    avg = pandas.concat(fracsurvives).groupby(cols[ : -1])
    if avgtype == 'mean':
        avg = avg.mean()
    elif avgtype == 'median':
        avg = avg.median()
    else:
        raise ValueError("invalid avgtype {0}".format(avgtype))
    return (avg.reset_index()
               .sort_values('mutfracsurvive', ascending=False)
               [cols]
               .reset_index(drop=True)
               )


if __name__ == '__main__':
    import doctest
    doctest.testmod()
