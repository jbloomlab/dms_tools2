"""
===========
diffsel
===========

Performs operations related to estimating differential selection.
"""

import math
import os
import io
import tempfile

import natsort
import numpy
import pandas
import scipy.stats
from dms_tools2 import CODONS, CODON_TO_AA


def beta_diversity(tidy_df, *, samplecol, sitecol, valcol,
                   index, dropdups=True):
    """Computes beta diversity of selection across sites.

    Args:
        `tidy_df` (pandas DataFrame)
            Tidy data frame with selection values.
        `samplecol` (str)
            Column in `tidy_df` giving different samples.
        `sitecol` (str)
            Column in `tidy_df` giving sites.
        `valcol` (str)
            Column in `tidy_df` giving value (e.g., 'positive_diffsel')
        `index` {'shannon', 'simpson'}
            Type of diversity index to use.
        `dropdups` (bool)
            If duplicate columns for a site, drop them? Otherwise error.

    Returns:
        For each sample, the diversity is calculated by examining the
        total fraction of selection (as defined in `valcol`) allocated
        to each site (as defined in `sitecol`) and considering that to
        represent the weighted "species" abundance. Alpha and gamma diversity
        are then calculated using the indicated index as in Jost,
        "Partitioning diversity into independent alpha and beta components"
        (https://pdfs.semanticscholar.org/88c4/67e65b42eff088d19354632f41df6ca5804d.pdf).
        The beta diversity is the ratio of gamma to alpha diversity:
        larger values indicate the samples are more diverse in sites of selection.

    As an example, first create a data frame for three samples.
    Sample `a` is similar to `b`, but both are quite different from `c`:

    >>> a_df = pandas.DataFrame({
    ...     'sample': 'a',
    ...     'site': [1, 2, 3, 4],
    ...     'positive_diffsel': [0.2, 5.1, 0.2, 0.3]})
    >>> b_df = pandas.DataFrame({
    ...     'sample': 'b',
    ...     'site': [1, 2, 3, 4],
    ...     'positive_diffsel': [0.3, 4.2, 0.2, 0.3]})
    >>> c_df = pandas.DataFrame({
    ...     'sample': 'c',
    ...     'site': [1, 2, 3, 4],
    ...     'positive_diffsel': [0.3, 0.2, 0.2, 4.3]})
    >>> tidy_df = pandas.concat([a_df, b_df, c_df], ignore_index=True)
    >>> tidy_df
       sample  site  positive_diffsel
    0       a     1               0.2
    1       a     2               5.1
    2       a     3               0.2
    3       a     4               0.3
    4       b     1               0.3
    5       b     2               4.2
    6       b     3               0.2
    7       b     4               0.3
    8       c     1               0.3
    9       c     2               0.2
    10      c     3               0.2
    11      c     4               4.3

    Calculate beta diversity of all three samples using Shannon index.
    The value is relatively high because sample `c` is a lot different
    than `a` and `b`:

    >>> numpy.allclose(beta_diversity(tidy_df,
    ...                      samplecol='sample',
    ...                      sitecol='site',
    ...                      valcol='positive_diffsel',
    ...                      index='shannon'),
    ...                1.4914, atol=1e-4)
    True

    If we repeat the same using the Simpson index we get a higher
    diversity because the large-selection sites (which are up-weighted by
    Simpson index relative to Shannon) are diferent between `c` and `a` / `b`:

    >>> numpy.allclose(beta_diversity(tidy_df,
    ...                      samplecol='sample',
    ...                      sitecol='site',
    ...                      valcol='positive_diffsel',
    ...                      index='simpson'),
    ...                1.6478, atol=1e-4)
    True

    If we calculate the beta diversity of just samples `a` and `b`, we
    get a smaller value because those samples are quite similar:

    >>> numpy.allclose(beta_diversity(tidy_df.query('sample in ["a", "b"]'),
    ...                      samplecol='sample',
    ...                      sitecol='site',
    ...                      valcol='positive_diffsel',
    ...                      index='shannon'),
    ...                1.0022, atol=1e-4)
    True

    Here the Simpson index gives lower diversity than the Shannon since
    `a` and `b` share the largest selection site:

    >>> numpy.allclose(beta_diversity(tidy_df.query('sample in ["a", "b"]'),
    ...                      samplecol='sample',
    ...                      sitecol='site',
    ...                      valcol='positive_diffsel',
    ...                      index='simpson'),
    ...                1.0008, atol=1e-4)
    True

    And of course the beta diversity of two identical samples is one:

    >>> numpy.allclose(beta_diversity(
    ...                      pandas.concat([a_df, a_df.assign(sample='a2')]),
    ...                      samplecol='sample',
    ...                      sitecol='site',
    ...                      valcol='positive_diffsel',
    ...                      index='shannon'),
    ...                1, atol=1e-4)
    True
    >>> numpy.allclose(beta_diversity(
    ...                      pandas.concat([a_df, a_df.assign(sample='a2')]),
    ...                      samplecol='sample',
    ...                      sitecol='site',
    ...                      valcol='positive_diffsel',
    ...                      index='simpson'),
    ...                1, atol=1e-4)
    True

    """
    cols = [samplecol, sitecol, valcol]
    if set(cols) > set(tidy_df.columns):
        raise ValueError(f"`tidy_df` lacks required columns: {cols}")
    if len(tidy_df[samplecol].unique()) < 2:
        raise ValueError('`tidy_df` must have at least two unique samples')
    df = (tidy_df
          [cols]
          .drop_duplicates()
          )
    if (len(df) != len(tidy_df)) and not dropdups:
        raise ValueError(f"`tidy_df` has duplicate entries in columns: {cols}")
    sites = set(df[sitecol].unique())
    if any(set(idf[sitecol]) != sites for _, idf in df.groupby(samplecol)):
        raise ValueError('not same set of sites for all samples in `tidy_df`')
    if any(len(sites) != len(idf[sitecol]) for _, idf in df.groupby(samplecol)):
        raise ValueError('duplicate sites for a sample in `tidy_df`.')

    if 'normval' in cols:
        raise ValueError('`tidy_df` cannot have column named "normval"')
    df = (df
          .assign(normval=lambda x: (x[valcol] /
                                     (x.groupby(samplecol)
                                      [valcol]
                                      .transform('sum')
                                      )
                                     )
                  )
          [[samplecol, sitecol, 'normval']]
          )

    if index == 'shannon':
        gamma_diversity = math.exp(scipy.stats.entropy(
            df
            .groupby(sitecol)
            .aggregate({'normval': 'mean'})
            ['normval']
            ))
        alpha_diversity = math.exp(
            df
            .groupby(samplecol)
            ['normval']
            .apply(scipy.stats.entropy)
            .reset_index()
            ['normval']
            .mean()
            )
    elif index == 'simpson':
        gamma_diversity = 1 / (
            df
            .groupby(sitecol)
            .aggregate({'normval': 'mean'})
            .assign(normval=lambda x: x['normval']**2)
            ['normval']
            .sum()
            )
        alpha_diversity = 1 / (
            df
            .assign(normval=lambda x: x['normval']**2)
            .groupby(samplecol)
            ['normval']
            .aggregate('sum')
            .reset_index()
            ['normval']
            .mean()
            )
    else:
        raise ValueError(f"invalid `index` of {index}")

    return gamma_diversity / alpha_diversity


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
            natsort.index_realsorted(tidy_df.site)))

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
        wtmask = (m['mutation'] == m['wildtype'])
        assert all(m[wtmask]['epsilon'] > 0), \
                "err counts of 0 for wildtype"
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


def df_read_filecols(df, filecols, *, order_sites=True):
    """Merges data frame with entries read from files.

    Designed to expand data frame listing CSV files with site or
    mutation-level selection information into a data frame listing
    the information in these CSV files.

    Args:
        `df` (pandas DataFrame)
            Each row gives files and associated information.
        `filecols` (list)
            List of columns in `df` that give filenames of CSV files
            to add to data frame. These CSV files cannot have column
            names already in `df`.
        `order_sites` (bool)
            Expect a `site` column, make it naturally sorted
            categorical variable, and add `isite` column
            that numbers sites 0, 1, ...

    Returns:
        A data frame where the entries in the files are now
        read as columns.

    >>> tf = tempfile.NamedTemporaryFile
    >>> with tf(mode='w') as sitediffsel1, tf(mode='w') as sitediffsel2, \\
    ...      tf(mode='w') as mutdiffsel1,  tf(mode='w') as mutdiffsel2:
    ...
    ...     # first sitediffsel file
    ...     _ = sitediffsel1.write('site,sitediffsel\\n'
    ...                            '1,3.2\\n'
    ...                            '-1,2.3\\n'
    ...                            '(HA2)1,0.1')
    ...     sitediffsel1.flush()
    ...
    ...     # first mutdiffsel file
    ...     _ = mutdiffsel1.write('site,wildtype,mutation,mutdiffsel\\n'
    ...                            '-1,A,C,-0.7\\n'
    ...                            '-1,A,G,3.0\\n'
    ...                            '1,C,A,1.2\\n'
    ...                            '1,C,G,2.0\\n'
    ...                            '(HA2)1,C,A,0.0\\n'
    ...                            '(HA2)1,C,G,0.1')
    ...     mutdiffsel1.flush()
    ...
    ...     # second sitediffsel file
    ...     _ = sitediffsel2.write('site,sitediffsel\\n'
    ...                            '(HA2)1,9.1\\n'
    ...                            '1,1.2\\n'
    ...                            '-1,0.3\\n')
    ...     sitediffsel2.flush()
    ...
    ...     # second mutdiffsel file
    ...     _ = mutdiffsel2.write('site,wildtype,mutation,mutdiffsel\\n'
    ...                            '-1,A,C,-0.2\\n'
    ...                            '-1,A,G,0.5\\n'
    ...                            '1,C,A,1.1\\n'
    ...                            '1,C,G,0.1\\n'
    ...                            '(HA2)1,C,A,9.0\\n'
    ...                            '(HA2)1,C,G,0.1')
    ...     mutdiffsel2.flush()
    ...
    ...     # data frame with files as columns
    ...     df = pandas.DataFrame({
    ...          'name':['sample_1', 'sample_2'],
    ...          'serum':['serum_1', 'serum_1'],
    ...          'sitediffsel_file':[sitediffsel1.name, sitediffsel2.name],
    ...          'mutdiffsel_file':[mutdiffsel1.name, mutdiffsel2.name]
    ...          })
    ...
    ...     # call df_read_filecols
    ...     (df_read_filecols(df, ['sitediffsel_file', 'mutdiffsel_file'])
    ...      .drop(columns=['sitediffsel_file', 'mutdiffsel_file']))
            name    serum    site  sitediffsel wildtype mutation  mutdiffsel  isite
    0   sample_1  serum_1       1          3.2        C        A         1.2      1
    1   sample_1  serum_1       1          3.2        C        G         2.0      1
    2   sample_1  serum_1      -1          2.3        A        C        -0.7      0
    3   sample_1  serum_1      -1          2.3        A        G         3.0      0
    4   sample_1  serum_1  (HA2)1          0.1        C        A         0.0      2
    5   sample_1  serum_1  (HA2)1          0.1        C        G         0.1      2
    6   sample_2  serum_1  (HA2)1          9.1        C        A         9.0      2
    7   sample_2  serum_1  (HA2)1          9.1        C        G         0.1      2
    8   sample_2  serum_1       1          1.2        C        A         1.1      1
    9   sample_2  serum_1       1          1.2        C        G         0.1      1
    10  sample_2  serum_1      -1          0.3        A        C        -0.2      0
    11  sample_2  serum_1      -1          0.3        A        G         0.5      0
    """
    if not len(df):
        raise ValueError('`df` has no rows')

    df_cols = set(df.columns)
    if 'dummy' in df_cols:
        raise ValueError('`df` has column named "dummy"')

    if not (set(filecols) <= df_cols):
        raise ValueError('`df` does not have all the `filecol` columns')

    df_filecols = []
    for row in df.iterrows():
        # get data frame of just row, with a dummy column for merging
        row_df = row[1].to_frame().transpose().assign(dummy=1)
        for col in filecols:
            filename = row_df.at[row[0], col]
            file_df = pandas.read_csv(filename).assign(dummy=1)
            if order_sites and 'site' not in file_df.columns:
                raise ValueError(f"no `site` column in {filename}")
            sharedcols = set(file_df.columns).intersection(df_cols)
            if sharedcols:
                raise ValueError(f"`df` and {filename} share columns "
                                 f"{sharedcols}")
            row_df = row_df.merge(file_df)
        df_filecols.append(row_df)

    df_filecols = (pandas.concat(df_filecols, ignore_index=True)
                   .drop('dummy', axis='columns')
                   )

    if order_sites:
        sites = natsort.realsorted(df_filecols['site'].unique())
        df_filecols = (
            df_filecols
            .assign(site=lambda x: pandas.Categorical(x['site'],
                                                      sites,
                                                      ordered=True),
                    isite=lambda x: x['site'].cat.codes)
            )

    return df_filecols



if __name__ == '__main__':
    import doctest
    doctest.testmod()
