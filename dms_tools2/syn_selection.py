"""
==============
syn_selection
==============

Identifies enrichment of synonymous codons.

"""

import pandas as pd
import numpy as np
import scipy.stats

from dms_tools2 import CODON_TO_AA

def syn_selection_by_codon(counts_pre, counts_post, pseudocount=0.5):
    """Identify sites with selection on synonymous codons.

    Runs two-tailed Fisher Exact test, which returns:
     - P-value reflecting the significance of codon_x enrichment.

    After calculating the Fisher P-value, a pseudocount is added
    for calculation of the odds ratio, which reflects the enrichment of
    codon_x after selection (relative to other synonymous codons).
     - The pseudocount removes 0's to avoid returning NA or inf odds ratios.
     - Rows where only one codon is represented pre-selection are dropped.

    Args:
        `counts_pre` (str or pandas.DataFrame)
            CSV file giving pre-selection codon counts with columns
            named 'site', 'wildtype', and list of codons. Can also
            be a pandas DataFrame containing the CSV file.
        `counts_post` (str or pandas.DataFrame)
            Like `counts_pre` but for the post-selection counts.
            CSV file giving post-selection codon counts in same format
            as `counts_pre`.
        'pseudocount' (float or int, default 0.5)
            Number to add to each codon count before calculating the odds
            ratio.

    Returns:
        A pandas DataFrame with the following columns:
         - 'site'
         - 'wildtype' : wildtype codon at site
         - 'codon' : codon we are analyzing at site
         - 'aa' : amino acid
         - 'codon_pre' : counts for codon of interest pre-selection
         - 'aa_pre' : counts for all codons for amino acid pre-selection
         - 'codon_post' : counts for codon of interest post-selection
         - 'aa_post' : counts for all codons for amino acid post-selection
         - 'odds_ratio' : enrichment of codon post-selection
         - 'P' : P-value calculated using Fisher's exact test

    Example:

    >>> pd.set_option('display.max_columns', None)  # display all columns
    >>> pd.set_option('expand_frame_repr', False)  # do not break lines
    >>> counts_pre = pd.DataFrame.from_records(
    ...         [(1, 'ATC', 5, 100, 10),
    ...          (2, 'ATT', 50, 10, 10),
    ...          ],
    ...         columns=['site', 'wildtype', 'ATT', 'ATC', 'ATA'],
    ...         )
    >>> counts_post = pd.DataFrame.from_records(
    ...         [(1, 'ATC', 5, 50, 75),
    ...          (2, 'ATT', 50, 9, 11),
    ...          ],
    ...         columns=['site', 'wildtype', 'ATT', 'ATC', 'ATA'],
    ...         )
    >>> syn_selection_by_codon(counts_pre, counts_post, 1)
       site wildtype codon aa  codon_pre  aa_pre  codon_post  aa_post  odds_ratio             P
    0     1      ATC   ATT  I          6     118           6      133    0.881890  1.000000e+00
    1     1      ATC   ATC  I        101     118          51      133    0.104685  1.798192e-15
    2     1      ATC   ATA  I         11     118          76      133   12.969697  7.500709e-17
    3     2      ATT   ATT  I         51      73          51       73    1.000000  1.000000e+00
    4     2      ATT   ATC  I         11      73          10       73    0.894661  1.000000e+00
    5     2      ATT   ATA  I         11      73          12       73    1.108793  1.000000e+00

    """

    # translate codons and calculate aa_count for pre- and post-selection datasets
    counts_df_list = []
    for df, df_type in [(counts_pre, 'pre'), (counts_post, 'post')]:
        df = (
         df.melt(id_vars=['site', 'wildtype'],
                 var_name='codon',
                 value_name=f"codon_{df_type}")
         .assign(aa=lambda x: x['codon'].map(CODON_TO_AA))  # translate codons
         [['site', 'wildtype', 'codon', 'aa', f"codon_{df_type}"]]
         )

        # group synonymous codons at each site
        aaGroups = df.groupby(['site', 'aa'])

        df = (
         df.assign(aa_count=aaGroups[f"codon_{df_type}"].transform(np.sum))
         # use in calculating significance of codon vs other synonymous codons
         .assign(syn_codons=lambda x: x['aa_count'] - x[f"codon_{df_type}"])
         .rename(columns={'aa_count': f"aa_{df_type}",
                          'syn_codons': f"syn_codons_{df_type}"})
         )

        counts_df_list.append(df)

    # merge counts_pre and counts_post
    df_merge = pd.merge(counts_df_list[0],
                        counts_df_list[1],
                        on=['site', 'wildtype', 'codon', 'aa'])

    # apply Fisher's exact test to codon x vs other synonymous codons
    df_merge['fishers'] = (
        df_merge.apply(lambda x: scipy.stats.fisher_exact(
            [[x.syn_codons_pre, x.syn_codons_post],
             [x.codon_pre, x.codon_post]]), axis=1)
    )

    # split fishers results into 'odds_ratio' and 'P', and drop 'odds_ratio'
    new_col_list = ['odds_ratio', 'P']
    for n, col in enumerate(new_col_list):
        df_merge[col] = df_merge['fishers'].map(lambda x: x[n])
    # drop columns for recalculation with pseudocount
    df_merge = df_merge.drop(['odds_ratio', 'fishers', 'aa_pre', 'syn_codons_pre',
                              'aa_post', 'syn_codons_post'], axis=1)

    for df_type in ['pre', 'post']:
        # add pseudocount
        df_merge[f'codon_{df_type}'] += pseudocount
        # recalculate aa and syn_codon sums with pseudocount
        aaGroups = df_merge.groupby(['site', 'aa'])
        df_merge = (
            df_merge.assign(aa_count=aaGroups[f"codon_{df_type}"].transform(np.sum))
                # use in calculating significance of codon vs other synonymous codons
                .assign(syn_codons=lambda x: x['aa_count'] - x[f"codon_{df_type}"])
                .rename(columns={'aa_count': f"aa_{df_type}",
                                 'syn_codons': f"syn_codons_{df_type}"})
        )

    # drop rows where syn_codons_pre = 0
    df_merge = df_merge[df_merge.syn_codons_pre != 0]

    # calculate odds ratio with pseudocount
    df_merge['odds_ratio'] = (
        df_merge.apply(lambda x:
                       ((x.codon_post / x.syn_codons_post) /
                        (x.codon_pre / x.syn_codons_pre)), axis=1)
    )

    # drop extra columns
    df = (df_merge
        .drop(['syn_codons_pre', 'syn_codons_post'], axis=1)
        .sort_values(['site', 'aa'])
        .reset_index(drop=True)
    [['site', 'wildtype', 'codon', 'aa', 'codon_pre', 'aa_pre', 'codon_post',
      'aa_post', 'odds_ratio', 'P']]
        )

    return df


if __name__ == '__main__':
    import doctest
    doctest.testmod()
