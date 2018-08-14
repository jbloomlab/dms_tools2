"""
===============
barcodes
===============

Operations for sequence barcodes and UMIs.
"""


import collections


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


if __name__ == '__main__':
    import doctest
    doctest.testmod()
