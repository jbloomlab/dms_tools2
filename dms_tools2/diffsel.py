"""
===========
diffsel
===========

Performs operations related to estimating differential selection.
"""

import os
import numpy
import pandas
from dms_tools2 import CODONS, CODON_TO_AA


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
        m['n{0}'.format(name)] = m_name['n']
        m['N{0}'.format(name)] = m_name['N']

    # error correction 
    if err is not None:
        m['epsilon'] = m['nerr'] / m['Nerr']
        assert all(m[m['mutation'] == m['wildtype']] > 0), \
                "err counts of 0 for wildtype"
        for name in ['sel', 'mock']:
            ncol = 'n{0}'.format(name)
            Ncol = 'N{0}'.format(name)
            m[ncol] = numpy.maximum(
                    (m[Ncol] * (m[ncol] / m[Ncol] - m['epsilon'])).where(
                            m['mutation'] == m['wildtype'],
                            m[ncol] / m['epsilon'])
                    )
            m[Ncol] = m.groupby('site').transform('sum')[ncol]
   
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



if __name__ == '__main__':
    import doctest
    doctest.testmod()
