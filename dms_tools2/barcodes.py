"""
===============
barcodes
===============

Operations for sequence barcodes and UMIs.
"""


import collections


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
        df = df.assign(sample=sample_col)
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



if __name__ == '__main__':
    import doctest
    doctest.testmod()
