"""
===========
diffsel
===========

Performs operations related to estimating differential selection.
"""

import os
import io
import tempfile

import natsort
import numpy
import pandas
from dms_tools2 import CODONS, CODON_TO_AA


def tidyToWide(tidy_df, valuecol):
    """Converts tidy `diffsel` data frame to wide form.

    The `diffsel` data frames returned by ``dms2_diffsel`` (and
    also other dataframes, such as the `fracsurvive` ones
    from ``dms_fracsurvive`` are in tidy form. This function
    converts them to wide form.

    Args:
        `tidy_df` (pandas DataFrame)
            Data frame in tidy form. Should have columns named
            `site`, `wildtype`, `mutation`, and something
            with the name matching `valuecol`.
        `valuecol` (string)
            Name of value column in `df`, such `diffsel` or
            `fracsurvive`.

    Returns:
        Wide form dataframe. Will have columns `site` (as string), 
        `wildtype`, and all characters (e.g., amino acids)
        for which values are given. Natural sorted by `site`.

    >>> tidy_df = pandas.read_csv(io.StringIO(
    ...     '''site wildtype mutation diffsel
    ...           3        A        D    -1.5 
    ...           3        A        C    10.1
    ...           2        A        C    10.1
    ...           1        C        D     9.5
    ...           1        C        A     0.2
    ...           2        A        D    -1.5'''),
    ...     delim_whitespace=True, index_col=False)
    >>> wide_df = tidyToWide(tidy_df, valuecol='diffsel')
    >>> print(wide_df.to_string(float_format=lambda x: '{0:.1f}'.format(x)))
      site   A    C    D wildtype
    0    1 0.2  0.0  9.5        C
    1    2 0.0 10.1 -1.5        A
    2    3 0.0 10.1 -1.5        A
    """
    assert isinstance(tidy_df, pandas.DataFrame)
    cols = ['site', 'wildtype', 'mutation', valuecol]
    assert set(cols) == set(tidy_df.columns), ('expected columns '
            '{0}\nactual columns {1}'.format(cols, tidy_df.columns))

    # make site a string
    tidy_df['site'] = tidy_df['site'].astype(str)

    # sort on site as here: https://stackoverflow.com/a/29582718
    tidy_df = tidy_df.reindex(index=natsort.order_by_index(tidy_df.index,
            natsort.index_natsorted(tidy_df.site, signed=True)))

    # convert to wide form, keeping wildtype identities
    tidy_df = tidy_df.set_index('site', drop=True)
    wt = tidy_df['wildtype']
    wide_df = (tidy_df.pivot(columns='mutation', values=valuecol)
                      .fillna(0.0)
                      .join(wt)
                      .reset_index()
                      )
    wide_df = wide_df.drop_duplicates().reset_index(drop=True)

    return wide_df


def computeMutDiffSel(sel, mock, countcharacters, pseudocount, 
        translate_to_aa, err=None, mincount=0):
    """Compute mutation differential selection.

    Args:
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
            estimating diffsel for amino acids, `False` otherwise.
        `err` (pandas.DataFrame or `None`)
            Optional error-control counts, in same format as `sel`.
        `mincount` (int >= 0)
            Report as `NaN` the mutdiffsel for any mutation in which
            neither `sel` nor `mock` has at least this many counts.

    Returns:
        A `pandas.DataFrame` with the mutation differential selection.
        Columns are `site`, `wildtype`, `mutation`, `mutdiffsel`.
    """
    assert pseudocount > 0

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
    m['nselP'] = m['nsel'] + pseudocount * numpy.maximum(1,
            m['Nsel'] / m['Nmock'])
    m['nmockP'] = m['nmock'] + pseudocount * numpy.maximum(1,
            m['Nmock'] / m['Nsel'])

    # create columns with wildtype counts
    m = m.merge(
            m.query('mutation==wildtype')[['site', 'nselP', 'nmockP']],
            on='site', suffixes=('', 'wt'))

    # compute mutdiffsel
    m['enrichment'] = (
            (m['nselP'] / m['nselPwt']) / (m['nmockP'] / m['nmockPwt'])
            ).where((m['mutation'] != m['wildtype']) &
                    ((m['nsel'] >= mincount) | (m['nmock'] >= mincount)),
                    numpy.nan)
    m['mutdiffsel'] = numpy.log2(m['enrichment'])

    return m[['site', 'wildtype', 'mutation', 'mutdiffsel']]


def mutToSiteDiffSel(mutdiffsel):
    """Computes sitediffsel from mutdiffsel.

    Args:
        `mutdiffsel` (pandas.DataFrame)
            Dataframe with mutdiffsel as returned by `computeMutDiffSel`

    Returns:
        The dataframe `sitediffsel`, which has the following columns:
            - `site`
            - `abs_diffsel`: sum of absolute values of mutdiffsel at site
            - `positive_diffsel`: sum of positive mutdiffsel at site
            - `negative_diffsel`: sum of negative mutdiffsel at site
            - `max_diffsel`: maximum mutdiffsel at site
            - `min_diffsel`: minimum mutdiffsel at site

    >>> mutdiffsel = (pandas.DataFrame({
    ...         'site':[1, 2, 3, 4],
    ...         'wildtype':['A', 'C', 'C', 'A'],
    ...         'A':[numpy.nan, 4.1, -0.1, numpy.nan],
    ...         'C':[-0.2, numpy.nan, numpy.nan, 0.3],
    ...         'G':[3.2, 0.1, -0.2, 0.1],
    ...         'T':[-0.2, 0.0, -0.1, 0.4],
    ...         })
    ...         .melt(id_vars=['site', 'wildtype'],
    ...               var_name='mutation', value_name='mutdiffsel')
    ...         .reset_index(drop=True)
    ...         )
    >>> sitediffsel = mutToSiteDiffSel(mutdiffsel)
    >>> all(sitediffsel.columns == ['site', 'abs_diffsel', 
    ...         'positive_diffsel', 'negative_diffsel', 
    ...         'max_diffsel', 'min_diffsel'])
    True
    >>> all(sitediffsel['site'] == [1, 2, 3, 4])
    True
    >>> numpy.allclose(sitediffsel['abs_diffsel'], [3.6, 4.2, 0.4, 0.8])
    True
    >>> numpy.allclose(sitediffsel['positive_diffsel'], [3.2, 4.2, 0, 0.8])
    True
    >>> numpy.allclose(sitediffsel['negative_diffsel'], [-0.4, 0, -0.4, 0])
    True
    >>> numpy.allclose(sitediffsel['max_diffsel'], [3.2, 4.1, 0, 0.4])
    True
    >>> numpy.allclose(sitediffsel['min_diffsel'], [-0.2, 0, -0.2, 0])
    True
    """
    aas = list(set(mutdiffsel.wildtype).union(mutdiffsel.mutation))
    sitediffsel = (mutdiffsel
                   .pivot(index='site', columns='mutation', values='mutdiffsel')
                   .fillna(0)
                   .assign(abs_diffsel=lambda x:
                                    x[aas].abs().sum(axis=1),
                           positive_diffsel=lambda x:
                                    x[aas][x[aas] >= 0].sum(axis=1),
                           negative_diffsel=lambda x:
                                    x[aas][x[aas] <= 0].sum(axis=1),
                           max_diffsel=lambda x:
                                    x[aas].max(axis=1),
                           min_diffsel=lambda x:
                                    x[aas].min(axis=1),
                          )
                   .reset_index()
                   [['site', 'abs_diffsel', 'positive_diffsel', 
                    'negative_diffsel', 'max_diffsel', 'min_diffsel']]
                   )
    return sitediffsel


def avgMutDiffSel(mutdiffselfiles, avgtype):
    """Gets mean or median mutation differential selection.

    Args:
        `mutdiffselfiles` (list)
            List of CSV files with mutdiffsel as returned by
            ``dms2_diffsel``.
        `avgtype` (str)
            Type of "average" to calculate. Possibilities:
                - `mean`
                - `median`

    Returns:
        A `pandas.DataFrame` containing the mean or median
        mutation differential selection.

    >>> tf = tempfile.NamedTemporaryFile
    >>> with tf(mode='w') as f1, tf(mode='w') as f2, tf(mode='w') as f3:
    ...     x = f1.write('site,wildtype,mutation,mutdiffsel\\n'
    ...                  '157,K,D,8.3\\n'
    ...                  '156,G,S,2.0')
    ...     f1.flush()
    ...     x = f2.write('site,wildtype,mutation,mutdiffsel\\n'
    ...                  '157,K,D,6.3\\n'
    ...                  '156,G,S,-1.1')
    ...     f2.flush()
    ...     x = f3.write('site,wildtype,mutation,mutdiffsel\\n'
    ...                  '157,K,D,2.2\\n'
    ...                  '156,G,S,0.0')
    ...     f3.flush()
    ...     mean = avgMutDiffSel([f1.name, f2.name, f3.name], 'mean')
    ...     median = avgMutDiffSel([f1.name, f2.name, f3.name], 'median')
    >>> (mean['site'] == [157, 156]).all()
    True
    >>> (median['site'] == [157, 156]).all()
    True
    >>> (mean['wildtype'] == ['K', 'G']).all()
    True
    >>> (median['wildtype'] == ['K', 'G']).all()
    True
    >>> (mean['mutation'] == ['D', 'S']).all()
    True
    >>> (median['mutation'] == ['D', 'S']).all()
    True
    >>> numpy.allclose(mean['mutdiffsel'], [5.6, 0.3])
    True
    >>> numpy.allclose(median['mutdiffsel'], [6.3, 0.0])
    True
    """
    diffsels = []
    cols = ['site', 'wildtype', 'mutation', 'mutdiffsel']
    for f in mutdiffselfiles:
        diffsels.append(pandas.read_csv(f)
                        [cols]
                        .sort_values(['site', 'wildtype', 'mutation'])
                        .reset_index()
                        )
        assert all([(diffsels[0][c] == diffsels[-1][c]).all() for c in
                ['site', 'wildtype', 'mutation']]), "files do not have same muts" 
    avg = pandas.concat(diffsels).groupby(cols[ : -1])
    if avgtype == 'mean':
        avg = avg.mean()
    elif avgtype == 'median':
        avg = avg.median()
    else:
        raise ValueError("invalid avgtype {0}".format(avgtype))
    return (avg.reset_index()
               .sort_values('mutdiffsel', ascending=False)
               [cols]
               .reset_index(drop=True)
               )


if __name__ == '__main__':
    import doctest
    doctest.testmod()
