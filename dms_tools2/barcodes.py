"""
===============
barcodes
===============

Operations for sequence barcodes and UMIs.
"""

import re
import collections
import itertools

import numpy
import pandas
import umi_tools.network


_umi_clusterer = umi_tools.network.UMIClusterer()


def almost_duplicated(barcodes, threshold=1):
    """Identifies nearly identical barcodes.

    This function minics the pandas `duplicated`
    function except it can also mark as duplicated almost
    but not exactly identical sequences.

    Args:
        `barcodes` (list or pandas Series)
            Lists all the barcodes, which should be strings
            of same length.
        `threshold` (int)
            Max number of mismatches for barcodes to be
            considered almost identical.

    Returns:
        A pandas Series of the same length as `barcodes`
        filled with bools that indicates whether a barcode
        is a near duplicate. For each group of barcodes within
        `threshold` of each other, only one barcode will have
        an entry of `False` (not duplicated) and the rest will
        be `True` (duplicated). The barcode listed as not duplicated
        is chosen as follows: (1) it has the most abundant barcode
        in the group, (2) among barcodes that are equivalent abundant
        we take the one listed first in `barcodes`.
        Groups are computed using the directional method of
        `umi_tools <https://github.com/CGATOxford/UMI-tools>`_.

    When `threshold` is zero, just like `pandas.duplicated`:

    >>> barcodes = ['CTC', 'ATG', 'ATG', 'ATA']
    >>> almost_duplicated(barcodes, threshold=0).equals(
    ...     pandas.Series([False, False, True, False]))
    True
    >>> pandas.Series(barcodes).duplicated().equals(
    ...     pandas.Series([False, False, True, False]))
    True

    But when `threshold` is 1, also identifies **almost** equal
    barcodes (in this case, ones that differ by <= 1 mutation:

    >>> almost_duplicated(barcodes, threshold=1).equals(
    ...     pandas.Series([False, False, True, True]))
    True

    All barcodes are within a distance of two, so all are
    near duplicates. Note how we mark as not duplicated
    (`False`) the first one listed among the most abundant
    barcode ('ATG'):

    >>> almost_duplicated(pandas.Series(barcodes), threshold=2).equals(
    ...     pandas.Series([True, False, True, True]))
    True

    """
    if threshold < 0:
        raise ValueError("`threshold` must be >= 0")
    if not isinstance(barcodes, pandas.Series):
        if isinstance(barcodes, collections.Iterable):
            barcodes = pandas.Series(barcodes)
        else:
            raise TypeError(f"`barcodes` invalid type {type(barcodes)}")
    barcodes = barcodes.str.encode('utf-8')

    counts = collections.Counter(barcodes)

    groups = []
    for group in _umi_clusterer(barcodes, counts, threshold):
        # keep only barcode(s) most abundant in group
        max_count = max(counts[barcode] for barcode in group)
        groups.append(frozenset(barcode for barcode in group
                      if counts[barcode] >= max_count))

    dups = []
    for barcode in barcodes.values:
        for g in groups:
            if barcode in g:
                dups.append(False)
                groups.remove(g)
                break
        else:
            dups.append(True)
    assert not groups

    return pandas.Series(dups, index=barcodes.index)


def fracIdentWithinBarcode(df, *, barcode_col='barcode',
        variant_col='variant', sample_col=None):
    """Gets fraction of identical variants within barcodes.

    This function is designed for the case when you have barcoded
    variants, and you want to determine how often the variants
    with the same barcode are identical. To calculate this, all
    variants with the same barcode are grouped, and for every pair
    of variants within a barcode we compute the fraction that 
    are identical (note that this fraction is the
    `Simpson diversity index without replacement <https://en.wikipedia.org/wiki/Diversity_index>`_.
    We then compute the average of this fraction over all pairs.

    In the case where we are sequencing barcoded genes, the square
    root of this fraction can be interpreted as the sequencing accuracy
    per variant.

    Args:
        `df` (pandas Data Frame)
            Holds barcodes, variants, and (optionally) sample names.
        `barcode_col` (str)
            Name of column holding barcodes.
        `variant_col` (str)
            Name of column holding variants (these are what we
            compare to see if they are identical).
        `sample_col` (str or None)
            Name of samples. We do the computation separately
            for each sample. Set to `None` if only one sample.

    Returns:
        A pandas data frame with columns named `fraction_identical`
        and `accuracy` that contain the computed fraction identical
        and accuracy for each sample.

    Here is an example. In this data frame, we see that
    of the four pairs (one for barcode *AT* and three for
    barcode *TA*), two are identical. As you can see
    the function correctly calculates that 50% are identical,
    and computes an accuracy that is the square root of 0.5:

    >>> df = pandas.DataFrame({
    ...     'barcode':['AT', 'AT', 'TG', 'TA', 'TA', 'TA'],
    ...     'variant':['v1', 'v1', 'v2', 'v3', 'v4', 'v3']})
    >>> fracIdentWithinBarcode(df)
       fraction_identical  accuracy
    0                 0.5  0.707107

    Now another example with two samples and non-standard column
    names. Note that we only get results for the two samples
    with barcodes found multiple times:

    >>> df = pandas.DataFrame({
    ...     'bc'    :['AT', 'AT', 'TG', 'TA', 'TA', 'TA', 'TA'],
    ...     'var'   :['v1', 'v1', 'v2', 'v3', 'v4', 'v3', 'v3'],
    ...     'sample':['s1', 's1', 's2', 's3', 's3', 's3', 's4']})
    >>> fracIdentWithinBarcode(df, sample_col='sample',
    ...     barcode_col='bc', variant_col='var')
      sample  fraction_identical  accuracy
    0     s1            1.000000   1.00000
    1     s3            0.333333   0.57735

    """
    if sample_col is None:
        sample_col = 'sample'
        drop_sample_col = True
        df = df.assign(sample='dummy')
    else:
        drop_sample_col = False

    for col in [barcode_col, variant_col, sample_col]:
        if col not in df.columns:
            raise ValueError(f"No columns {col} in df")

    result = (
        df
        # get just sequences that have a barcode found multiple times
        .assign(barcode_counts=1)
        .assign(barcode_counts=lambda x:
            x.groupby([sample_col, barcode_col]).transform('count'))
        .query('barcode_counts > 1')
        # within each barcode, count number of sequences of each mutation combo
        .assign(dummy=1)
        .groupby([sample_col, barcode_col, 'barcode_counts', variant_col])
        .dummy
        .count()
        .reset_index(name='sequence_counts')
        # compute Simpson diversity without replacement for each barcode
        .groupby([sample_col, barcode_col, 'barcode_counts'])
        .apply(lambda x: (x.sequence_counts * (x.sequence_counts - 1) /
                      x.barcode_counts / (x.barcode_counts - 1)).sum())
        .reset_index(name='simpson_diversity')
        # compute weighted average of fraction identical across all pairs
        .assign(npairs=lambda x: x.barcode_counts * (x.barcode_counts - 1) / 2,
            weighted_diversity=lambda x: x.npairs * x.simpson_diversity)
        .groupby(sample_col)
        .apply(lambda x: x.weighted_diversity.sum() / x.npairs.sum())
        .reset_index(name='fraction_identical')
        # estimate accuracy as square root of fraction identical
        .assign(accuracy=lambda x: numpy.sqrt(x.fraction_identical))
        )

    if drop_sample_col:
        result = result.drop(sample_col, axis='columns')

    return result


def simpleConsensus(df, *,
        barcode_col='barcode', substitution_col='substitutions',
        insertion_col='insertions', deletion_col='deletions',
        sample_col=None, max_sub_diffs=1, max_indel_diffs=2,
        max_minor_muts=1):
    """Simple method to get consensus of mutations within barcode.

    Args:
        `df` (pandas Data Frame)
            Holds variants and their barcodes. Each row gives a 
            sequence variant and its barcode. There need to be
            columns with the names given by the next four arguments
            described below.
        `barcode_col` (str)
            Name of column holding barcodes.
        `substitution_col` (str)
            Name of column holding substitutions as list of strings.
        `insertion_col` (str)
            Name of column holding insertions as list of strings.
        `deletion_col` (str)
            Name of column holding insertions as list of strings.
        `sample_col` (`None` or str)
            If we have multiple samples, analyze each barcode only
            within its sample. In that case, `sample_col` should be
            name of column giving sample name.
        `max_sub_diffs` (int)
            Drop any barcode where any variant differs from all other
            variants for that barcode by more than this many substitution
            (point mutation) differences.
        `max_indel_diffs` (int)
            Drop any barcode where any variant differs from all other
            variants for that barcode by more than this many indel
            (insertion or deletion) differences.
        `max_minor_muts` (int)
            Drop any barcode where there is a minor (non-consensus)
            mutation found more than this many times.

    Returns:
        The 2-tuple `(consensus, dropped)`. These are each data frames:

            - `consensus` is a new data frame with a row for each
              barcode for which we could call the consensus. The
              columns have the same names as `barcode_col`, 
              `substitution_col`, `insertion_col`, `deletion_col`,
              and (optionally) `sample_col`--but the the three
              columns for the mutations now just list the **consensus**
              mutations of that type. In addition, there is a new
              column called `nsequences` that gives the number of
              sequences supporting the call of that barcode.

            - `dropped` simply contains all rows in the original `df`
              that correspond to sequences that were dropped due to
              `max_diffs` or `max_minor_muts` not being satisfied.
              There is also a column called "drop_reason" that is
              "excess diffs" if a sequence is dropped due to the
              `max_diffs` filter (too many differences between
              variants) or "excess minor muts" if a sequence is
              dropped due to the `max_minor_muts` filter (high
              frequency non-consensus mutations).

    The approach is as follows:

      1. Group all variants within the sample barcode and sample.

      2. If there are multiple sequences, check if any of them differ
         from all the others by more than `max_diffs` mutations total
         (taking substitutions, insertions, and deletions together).
         If so, drop the entire barcode. The reason is that if there
         are variants that are very different, it becomes likely that
         it isn't just sequencing error, but rather something is wrong
         with that barcode (multiple variants with same barcode or
         strand exchange).

      3. Take the consensus of the sequences, which means keeping
         mutations that are present in **greater** than half of the
         variants. Note that this calling scheme means that the
         consensus being called is dependent on the reference used
         to call the mutations, which is an important caveat if you
         are calling variants relative to multiple different parent
         sequences.

      4. If there are any minor mutations (mutations not in consensus)
         that are present in more than `max_minor_muts` variants or missing
         from more than `max_minor_muts` variants, then
         drop that barcode. The reason is that recurring minor mutations
         also suggest some problem more complex than sequencing error
         that may render the whole barcode family invalid.

    Note that this method returns a consensus even if there is just
    one sequence for the barcode (in that case, this sequence is
    the consensus). This is fine--if you want to get consensus calls
    that are more strongly supported, simply filter the returned
    `consensus` data frame for lager values of `nsequences`, as the
    more sequences that support a barcode call the more accurate is
    is expected to be.

    Here is an example:

    >>> df = pandas.DataFrame([
    ...     ('s1', 'AG', ['A2C'], [], ['del5to7']),
    ...     ('s1', 'AG', ['A2C'], [], []),
    ...     ('s1', 'TA', ['G3A'], ['ins4len3'], []),
    ...     ('s2', 'TA', ['C5A', 'T6C'], [], []),
    ...     ('s2', 'TA', ['T6C'], ['ins5len1'], []),
    ...     ('s2', 'TA', ['T6C'], [], []),
    ...     ('s2', 'GG', [], [], ['del1to4']),
    ...     ('s2', 'GG', ['A1C'], [], []),
    ...     ('s2', 'AA', [], [], []),
    ...     ('s2', 'AA', [], [], []),
    ...     ('s2', 'AA', ['T6C'], [], []),
    ...     ('s2', 'AA', ['T6C'], [], [])],
    ...     columns=['sample', 'barcode', 'substitutions',
    ...              'insertions', 'deletions']
    ...     )
    >>> consensus, dropped = simpleConsensus(df, sample_col='sample')
    >>> consensus
      sample barcode substitutions  insertions deletions  nsequences
    0     s1      AG         [A2C]          []        []           2
    1     s1      TA         [G3A]  [ins4len3]        []           1
    2     s2      TA         [T6C]          []        []           3
    >>> pandas.set_option('display.max_columns', 10)
    >>> dropped
      sample barcode substitutions insertions  deletions        drop_reason
    0     s2      GG            []         []  [del1to4]       excess diffs
    1     s2      GG         [A1C]         []         []       excess diffs
    2     s2      AA            []         []         []  excess minor muts
    3     s2      AA            []         []         []  excess minor muts
    4     s2      AA         [T6C]         []         []  excess minor muts
    5     s2      AA         [T6C]         []         []  excess minor muts
    """
    if sample_col is None:
        sample_col = 'sample'
        df = df.assign(sample_col='dummy')
        drop_sample_col = True
    else:
        drop_sample_col = False

    mut_cols = [substitution_col, insertion_col, deletion_col]
    all_cols = [sample_col, barcode_col] + mut_cols

    if not all([col in df.columns for col in all_cols]):
        raise ValueError(f"Cannot find column {col}")


    # make sure no mutations duplicated, otherwise approach below fails
    for col in mut_cols:
        duplicated = df[col].apply(len) - df[col].apply(set).apply(len)
        if duplicated.any():
            raise ValueError(f"duplicated {col}:\n"
                             f"{df[col][duplicated > 0]}")

    dropped = []
    consensus = []

    for (sample, barcode), g in df[all_cols].reset_index(drop=True).groupby(
            [sample_col, barcode_col]):

        nseqs = len(g)

        if nseqs == 1:
            consensus.append(g.values[0].tolist() + [nseqs])
            continue

        consensus_failed = False
        # are max_sub_diffs and max_indel_diffs satisfied?
        for difftype, mut_cols, max_diffs in [
                ('substitutions', [substitution_col], max_sub_diffs),
                ('indels', [insertion_col, deletion_col], max_indel_diffs)]:
            min_variant_diffs = collections.defaultdict(lambda: max_diffs + 1)
            for v1, v2 in itertools.combinations(g.itertuples(), 2):
                i1 = getattr(v1, 'Index')
                i2 = getattr(v2, 'Index')
                n_diffs = sum([len(
                        set(getattr(v1, col)).symmetric_difference(
                        set(getattr(v2, col))))
                        for col in mut_cols])
                min_variant_diffs[i1] = min(
                        min_variant_diffs[i1], n_diffs)

            if nseqs > 1 and any(
                    [d > max_diffs for d in min_variant_diffs.values()]):
                # need to add to `dropped` because of max_diffs failing
                dropped.append(g.assign(drop_reason=f"excess {difftype}"))
                consensus_failed = True
                break
        if consensus_failed:
            continue

        # get consensus and see if `max_minor_muts` is satisfied
        g_consensus = [sample, barcode]
        for col in mut_cols:
            counts = collections.Counter(
                    itertools.chain.from_iterable(g[col]))
            if any([max_minor_muts < count < (nseqs - max_minor_muts)
                   for count in counts.values()]):
                consensus_failed = True
                break
            else:
                col_consensus = [mut for mut, c in counts.items()
                        if c > 0.5 * nseqs]
                # order mutations based on first number in string
                n_col_consensus = []
                for mut in col_consensus:
                    m = re.search('(\-{0,1}\d+)', mut)
                    if m is None:
                        n_col_consensus.append((math.nan, mut))
                    else:
                        n_col_consensus.append((int(m.group()), mut))
                g_consensus.append([mut for n, mut in
                        sorted(n_col_consensus)])
        if consensus_failed:
            # need to add to dropped
            dropped.append(g.assign(drop_reason="excess minor muts"))
        else:
            consensus.append(g_consensus + [nseqs])

    consensus = pandas.DataFrame(consensus,
            columns=all_cols + ['nsequences'])
    dropped = pandas.concat(dropped).sort_index().reset_index(drop=True)

    if drop_sample_col:
        dropped = dropped.drop(sample_col, axis='columns')
        consensus = consensus.drop(sample_col, axis='columns')

    return (consensus, dropped)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
