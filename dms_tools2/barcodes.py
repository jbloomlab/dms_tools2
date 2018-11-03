"""
===============
barcodes
===============

Operations for sequence barcodes and UMIs.
"""

import re
import os
import collections
import itertools
import tempfile

import numpy
import pandas
import regex
import umi_tools.network

# use plotnine for plotting
from plotnine import *

import dms_tools2.plot
from dms_tools2.plot import COLOR_BLIND_PALETTE_GRAY
import dms_tools2.utils
import dms_tools2.pacbio
from dms_tools2 import CODON_TO_AA

_umi_clusterer = umi_tools.network.UMIClusterer()


def almost_duplicated(barcodes, threshold=1):
    """Identifies nearly identical barcodes.

    This function mimics the pandas `duplicated`
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
        variant_col='variant', library_col=None):
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
            Holds barcodes, variants, and (optionally) library names.
        `barcode_col` (str)
            Name of column holding barcodes.
        `variant_col` (str)
            Name of column holding variants (these are what we
            compare to see if they are identical).
        `library_col` (str or None)
            Name of libary. We do the computation separately
            for each library. Set to `None` if only one library.

    Returns:
        A pandas data frame with columns named `fraction_identical`
        and `accuracy` that contain the computed fraction identical
        and accuracy for each library.

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

    Now another example with two libraries and non-standard column
    names. Note that we only get results for the two libraries
    with barcodes found multiple times:

    >>> df = pandas.DataFrame({
    ...     'bc'    :['AT', 'AT', 'TG', 'TA', 'TA', 'TA', 'TA'],
    ...     'var'   :['v1', 'v1', 'v2', 'v3', 'v4', 'v3', 'v3'],
    ...     'library':['s1', 's1', 's2', 's3', 's3', 's3', 's4']})
    >>> fracIdentWithinBarcode(df, library_col='library',
    ...     barcode_col='bc', variant_col='var')
      library  fraction_identical  accuracy
    0      s1            1.000000   1.00000
    1      s3            0.333333   0.57735

    """
    if library_col is None:
        library_col = 'library'
        drop_library_col = True
        df = df.assign(library='dummy')
    else:
        drop_library_col = False

    for col in [barcode_col, variant_col, library_col]:
        if col not in df.columns:
            raise ValueError(f"No columns {col} in df")

    result = (
        df
        # get just sequences that have a barcode found multiple times
        .assign(barcode_counts=1)
        .assign(barcode_counts=lambda x:
            x.groupby([library_col, barcode_col]).transform('count'))
        .query('barcode_counts > 1')
        # within each barcode, count number of sequences of each mutation combo
        .assign(dummy=1)
        .groupby([library_col, barcode_col, 'barcode_counts', variant_col])
        .dummy
        .count()
        .reset_index(name='sequence_counts')
        # compute Simpson diversity without replacement for each barcode
        .groupby([library_col, barcode_col, 'barcode_counts'])
        .apply(lambda x: (x.sequence_counts * (x.sequence_counts - 1) /
                      x.barcode_counts / (x.barcode_counts - 1)).sum())
        .reset_index(name='simpson_diversity')
        # compute weighted average of fraction identical across all pairs
        .assign(npairs=lambda x: x.barcode_counts * (x.barcode_counts - 1) / 2,
            weighted_diversity=lambda x: x.npairs * x.simpson_diversity)
        .groupby(library_col)
        .apply(lambda x: x.weighted_diversity.sum() / x.npairs.sum())
        .reset_index(name='fraction_identical')
        # estimate accuracy as square root of fraction identical
        .assign(accuracy=lambda x: numpy.sqrt(x.fraction_identical))
        )

    if drop_library_col:
        result = result.drop(library_col, axis='columns')

    return result


def simpleConsensus(df, *,
        barcode_col='barcode', substitution_col='substitutions',
        insertion_col='insertions', deletion_col='deletions',
        library_col=None, max_sub_diffs=1, max_indel_diffs=2,
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
        `library_col` (`None` or str)
            If we have multiple libraries, analyze each barcode only
            within its library. In that case, `library_col` should be
            name of column giving library name.
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
              and (optionally) `library_col`--but the the three
              columns for the mutations now just list the **consensus**
              mutations of that type. In addition, there is a new
              column called "variant_call_support" that gives the number
              of sequences supporting the call of that barcode.

            - `dropped` simply contains all rows in the original `df`
              that correspond to sequences that were dropped due to
              `max_diffs` or `max_minor_muts` not being satisfied.
              There is also a column called "drop_reason" that gives
              the reason that the barcode was dropped.

    The approach is as follows:

      1. Group all variants by library and barcode.

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
    `consensus` data frame for larger values of `variant_call_support`,
    as the more sequences that support a call the more accurate it
    is expected to be.

    Here is an example:

    >>> df = pandas.DataFrame([
    ...     ('s1', 'AG', ['A2C'], [], ['del5to7']),
    ...     ('s1', 'AG', ['A2C'], [], []),
    ...     ('s1', 'TA', ['G3A'], ['ins4len3'], []),
    ...     ('s2', 'TA', ['C5A', 'T6C'], [], []),
    ...     ('s2', 'TA', ['T6C'], ['ins5len1'], []),
    ...     ('s2', 'TA', ['T6C'], [], []),
    ...     ('s2', 'TG', ['T6A'], [], []),
    ...     ('s2', 'TG', ['A2G'], [], []),
    ...     ('s2', 'GG', [], [], ['del1to4']),
    ...     ('s2', 'GG', ['A1C'], [], []),
    ...     ('s2', 'AA', [], [], []),
    ...     ('s2', 'AA', [], [], []),
    ...     ('s2', 'AA', ['T6C'], [], []),
    ...     ('s2', 'AA', ['T6C'], [], []),
    ...     ('s3', 'AA', ['T6G'], ['ins1len1'], ['del1to2']),
    ...     ('s3', 'AA', ['T6G'], [], ['del5to7'])],
    ...     columns=['library', 'barcode', 'substitutions',
    ...              'insertions', 'deletions']
    ...     )
    >>> consensus, dropped = simpleConsensus(df, library_col='library')
    >>> consensus
      library barcode substitutions  insertions deletions  variant_call_support
    0      s1      AG         [A2C]          []        []                     2
    1      s1      TA         [G3A]  [ins4len3]        []                     1
    2      s2      GG            []          []        []                     2
    3      s2      TA         [T6C]          []        []                     3
    >>> pandas.set_option('display.max_columns', 10)
    >>> pandas.set_option('display.width', 500)
    >>> dropped
      library barcode substitutions  insertions  deletions           drop_reason
    0      s2      TG         [T6A]          []         []  excess substitutions
    1      s2      TG         [A2G]          []         []  excess substitutions
    2      s2      AA            []          []         []     excess minor muts
    3      s2      AA            []          []         []     excess minor muts
    4      s2      AA         [T6C]          []         []     excess minor muts
    5      s2      AA         [T6C]          []         []     excess minor muts
    6      s3      AA         [T6G]  [ins1len1]  [del1to2]         excess indels
    7      s3      AA         [T6G]          []  [del5to7]         excess indels
    """
    if library_col is None:
        library_col = 'library'
        df = df.assign(library_col='dummy')
        drop_library_col = True
    else:
        drop_library_col = False

    mut_cols = [substitution_col, insertion_col, deletion_col]
    all_cols = [library_col, barcode_col] + mut_cols

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

    for (library, barcode), g in df[all_cols].reset_index(drop=True).groupby(
            [library_col, barcode_col]):

        nseqs = len(g)

        if nseqs == 1:
            consensus.append(g.values[0].tolist() + [nseqs])
            continue

        consensus_failed = False
        # are max_sub_diffs and max_indel_diffs satisfied?
        for difftype, diff_cols, max_diffs in [
                ('substitutions', [substitution_col], max_sub_diffs),
                ('indels', [insertion_col, deletion_col], max_indel_diffs)]:
            min_variant_diffs = collections.defaultdict(lambda: max_diffs + 1)
            for v1, v2 in itertools.combinations(g.itertuples(), 2):
                i1 = getattr(v1, 'Index')
                i2 = getattr(v2, 'Index')
                n_diffs = sum([len(
                        set(getattr(v1, col)).symmetric_difference(
                        set(getattr(v2, col))))
                        for col in diff_cols])
                min_variant_diffs[i1] = min(
                        min_variant_diffs[i1], n_diffs)

            if nseqs > 1 and any(
                    [d > max_diffs for d in min_variant_diffs.values()]):
                # need to add to `dropped` because of max_diffs failing
                dropped.append(
                        g.assign(drop_reason=f"excess {difftype}")
                        )
                consensus_failed = True
                break
        if consensus_failed:
            continue

        # get consensus and see if `max_minor_muts` is satisfied
        g_consensus = [library, barcode]
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
            dropped.append(
                    g.assign(drop_reason="excess minor muts")
                    )
        else:
            consensus.append(g_consensus + [nseqs])

    consensus = pandas.DataFrame(consensus,
            columns=all_cols + ['variant_call_support'])
    if dropped:
        dropped = pandas.concat(dropped).sort_index().reset_index(drop=True)
    else:
        dropped = pandas.DataFrame()

    if drop_library_col:
        dropped = dropped.drop(library_col, axis='columns')
        consensus = consensus.drop(library_col, axis='columns')

    return (consensus, dropped)


class IlluminaBarcodeParser:
    """Parser for Illumina barcodes.

    The barcodes should be read by R1 and optionally R2.
    The arrangement of elements is shown below::

        5'-[R2_start]-upstream-barcode-downstream-R1_start-3'

    R1 anneals downstream of the barcode and reads backwards. If 
    R2 is used, it anneals upstream of the barcode and reads forward.
    There can be sequences (`upstream` and `downstream`) on either
    side of the barcode: `downstream` must fully cover the
    region between where R1 starts and the barcode, and if you are
    using R2 then `upstream` must fully cover the region between
    where R2 starts and the barcode. However, it is fine if R1
    reads backwards past `upstream`, and if `R2` reads forward
    past `downstream`.

    Args:
        `bclen` (int)
            Length of the barcode.
        `upstream` (str)
            Sequence upstream of the barcode.
        `downstream` (str)
            Sequence downstream of barcode.
        `upstream_mismatch` (int)
            Max number of mismatches allowed in `upstream`.
        `downstream_mismatch` (int)
            Like `upstream_mismatches` but for `downstream`.
        `valid_barcodes` (`None` or iterable such as list, Series)
            If not `None`, only retain barcodes listed here.
            Use if you know the set of possible valid barcodes.
        `rc_barcode` (bool)
            Parse the reverse complement of the barcode (the
            orientation read by R1).
        `minq` (int)
            Require at least this quality score for all bases
            in barcode.
        `chastity_filter` (bool)
            Drop any reads that fail Illumina chastity filter.
        `list_all_valid_barcodes` (bool)
            If using `valid_barcodes`, then barcode sets returned
            by :class:`IlluminaBarcodeParser.parse` includes all
            valid barcodes even if no counts.

    To use, first initialize a :class:`IlluminaBarcodeParser`, then 
    parse barcodes using :class:`IlluminaBarcodeParser.parse`.
    Barcodes are retained as valid only if R1 and R2 agree at every
    nucleotide in barcode and if at each site at least one read has
    a quality of at least `minq`.

    Here is an example. Imagine we are parsing 4 nucleotide barcodes
    that have the following construction::

        5'-[R2 binding site]-ACATGA-NNNN-GACT-[R1 binding site]-3'

    First, we initialize an appropriate :class:`IlluminaBarcodeParser`:

    >>> parser = IlluminaBarcodeParser(
    ...              bclen=4,
    ...              upstream='ACATGA',
    ...              downstream='GACT'
    ...              )

    Now we write some test FASTQ files. We write valid test
    reads and some invalid reads. The header for each read
    explains why it is valid / invalid. We use quality scores
    of ``?`` (30) or ``+`` (10) for high- and low-quality bases:

    >>> r1file = '_temp_R1.fastq'
    >>> r2file = '_temp_R2.fastq'
    >>> with open(r1file, 'w') as f1, open(r2file, 'w') as f2:
    ...
    ...     # valid TACG barcode, full flanking regions
    ...     _ = f1.write(
    ...         '@valid_CGTA_barcode_full_flanking_region\\n'
    ...         'AGTCCGTATCATGT\\n'
    ...         '+\\n'
    ...         '??????????????\\n')
    ...     _ = f2.write(
    ...         '@valid_CGTA_barcode_full_flanking_region\\n'
    ...         'ACATGATACGGACT\\n'
    ...         '+\\n'
    ...         '??????????????\\n')
    ...
    ...     # valid CGTA barcode, partial flanking regions
    ...     _ = f1.write(
    ...         '@valid_CGTA_barcode_partial_flanking_region\\n'
    ...         'AGTCCGTATCAT\\n'
    ...         '+\\n'
    ...         '??????????????\\n')
    ...     _ = f2.write(
    ...         '@valid_CGTA_barcode_partial_flanking_region\\n'
    ...         'ACATGATACG\\n'
    ...         '+\\n'
    ...         '??????????????\\n')
    ...
    ...     # valid GCCG barcode, extended flanking regions
    ...     _ = f1.write(
    ...         '@valid_GCCG_barcode_extended_flanking_region\\n'
    ...         'AGTCGCCGTCATGTTAC\\n'
    ...         '+\\n'
    ...         '??????????????\\n')
    ...     _ = f2.write(
    ...         '@valid_GCCG_barcode_extended_flanking_region\\n'
    ...         'ACATGACGGCGACTGAC\\n'
    ...         '+\\n'
    ...         '??????????????\\n')
    ...
    ...     # AAGT barcode in R1 but R2 differs
    ...     _ = f1.write(
    ...         '@AAGT_R1_barcode_but_R2_differs\\n'
    ...         'AGTCAAGTTCATGT\\n'
    ...         '+\\n'
    ...         '??????????????\\n')
    ...     _ = f2.write(
    ...         '@AAGT_R1_barcode_but_R2_differs\\n'
    ...         'ACATGAACTAGACT\\n'
    ...         '+\\n'
    ...         '??????????????\\n')
    ...
    ...     # same site low quality in R1 and R2
    ...     _ = f1.write(
    ...         '@low_quality_site_in_R1_and_R2\\n'
    ...         'AGTCCGTATCATGT\\n'
    ...         '+\\n'
    ...         '?????+????????\\n')
    ...     _ = f2.write(
    ...         '@low_quality_site_in_R1_and_R2\\n'
    ...         'ACATGATACGGACT\\n'
    ...         '+\\n'
    ...         '????????+?????\\n')
    ...
    ...     # different site low quality in R1 and R2
    ...     _ = f1.write(
    ...         '@AGTA_with_low_quality_site_in_R1\\n'
    ...         'AGTCAGTATCATGT\\n'
    ...         '+\\n'
    ...         '?????+????????\\n')
    ...     _ = f2.write(
    ...         '@AGTA_with_low_quality_site_in_R1\\n'
    ...         'ACATGATACTGACT\\n'
    ...         '+\\n'
    ...         '?????????+????\\n')
    ...
    ...     # N in barcode
    ...     _ = f1.write(
    ...         '@N_in_barcode\\n'
    ...         'AGTCCGTNTCATGT\\n'
    ...         '+\\n'
    ...         '??????????????\\n')
    ...     _ = f2.write(
    ...         '@N_in_barcode\\n'
    ...         'ACATGATACGGACT\\n'
    ...         '+\\n'
    ...         '??????????????\\n')
    ...
    ...     # GGAG barcode, one mismatch in each flanking region
    ...     _ = f1.write(
    ...         '@GGAG_barcode_one_mismatch_per_flank\\n'
    ...         'GGTCGGAGTCATGA\\n'
    ...         '+\\n'
    ...         '??????????????\\n')
    ...     _ = f2.write(
    ...         '@GGAG_barcode_one_mismatch_per_flank\\n'
    ...         'TCATGACTCCGACG\\n'
    ...         '+\\n'
    ...         '??????????????\\n')
    ...
    ...     # GGAG barcode, two mismatch in a flanking region
    ...     _ = f1.write(
    ...         '@GGAG_barcode_two_mismatch_in_a_flank\\n'
    ...         'GGTCGGAGTCATAA\\n'
    ...         '+\\n'
    ...         '??????????????\\n')
    ...     _ = f2.write(
    ...         '@GGAG_barcode_two_mismatch_in_a_flank\\n'
    ...         'TCATGACTCCGACG\\n'
    ...         '+\\n'
    ...         '??????????????\\n')


    Now parse the barcodes using both the R1 and R2 files:

    >>> barcodes, fates = parser.parse(r1file, r2file)
    >>> print(barcodes.to_csv(sep=' ', index=False).strip())
    barcode count
    CGTA 2
    AGTA 1
    GCCG 1
    >>> print(fates.to_csv(sep=' ', index=False).strip())
    fate count
    "valid barcode" 4
    "unparseable barcode" 3
    "R1 / R2 disagree" 1
    "low quality barcode" 1

    Now we parse just using R1. We gain the barcode where R1 and
    R2 disagree, but lose the one where R1 is low quality at a
    position where R2 is OK:

    >>> barcodes, fates = parser.parse(r1file)
    >>> print(barcodes.to_csv(sep=' ', index=False).strip())
    barcode count
    CGTA 2
    AAGT 1
    GCCG 1
    >>> print(fates.to_csv(sep=' ', index=False).strip())
    fate count
    "valid barcode" 4
    "unparseable barcode" 3
    "low quality barcode" 2

    Now create a parser that allows a mismatch in each flanking
    region, and check that we recover a "GGAG" barcode:

    >>> parser_mismatch = IlluminaBarcodeParser(
    ...              bclen=4,
    ...              upstream='ACATGA',
    ...              downstream='GACT',
    ...              upstream_mismatch=1,
    ...              downstream_mismatch=1,
    ...              )
    >>> barcodes_mismatch, fates_mismatch = parser_mismatch.parse(r1file, r2file)
    >>> print(barcodes_mismatch.to_csv(sep=' ', index=False).strip())
    barcode count
    CGTA 2
    AGTA 1
    GCCG 1
    GGAG 1
    >>> print(fates_mismatch.to_csv(sep=' ', index=False).strip())
    fate count
    "valid barcode" 5
    "unparseable barcode" 2
    "R1 / R2 disagree" 1
    "low quality barcode" 1

    Now parse the barcodes using `valid_barcodes` to set a
    barcode whitelist:

    >>> parser_wl = IlluminaBarcodeParser(
    ...              bclen=4,
    ...              upstream='ACATGA',
    ...              downstream='GACT',
    ...              valid_barcodes={'CGTA', 'AGTA', 'TAAT'}
    ...              )
    >>> barcodes_wl, fates_wl = parser_wl.parse(r1file, r2file)
    >>> print(barcodes_wl.to_csv(sep=' ', index=False).strip())
    barcode count
    CGTA 2
    AGTA 1
    TAAT 0
    >>> print(fates_wl.to_csv(sep=' ', index=False).strip())
    fate count
    "unparseable barcode" 3
    "valid barcode" 3
    "R1 / R2 disagree" 1
    "invalid barcode" 1
    "low quality barcode" 1

    Remove the test FASTQ files:

    >>> os.remove(r1file)
    >>> os.remove(r2file)
    """

    #: valid nucleotide characters
    VALID_NTS = 'ACGTN'

    def __init__(self, *, bclen,
            upstream='', downstream='',
            upstream_mismatch=0, downstream_mismatch=0,
            valid_barcodes=None, rc_barcode=True, minq=20,
            chastity_filter=True, list_all_valid_barcodes=True):
        """See main class doc string."""

        # first make all arguments into attributes 
        self.bclen = bclen
        if re.match(f"^[{self.VALID_NTS}]*$", upstream):
            self.upstream = upstream
        else:
            raise ValueError(f"invalid chars in upstream {upstream}")
        if re.match(f"^[{self.VALID_NTS}]*$", downstream):
            self.downstream = downstream
        else:
            raise ValueError(f"invalid chars in downstream {downstream}")
        self.upstream_mismatch = upstream_mismatch
        self.downstream_mismatch = downstream_mismatch
        self.valid_barcodes = valid_barcodes
        if self.valid_barcodes:
            self.valid_barcodes = set(self.valid_barcodes)
        self.minq = minq
        self.rc_barcode = rc_barcode
        self.chastity_filter = chastity_filter
        self.list_all_valid_barcodes = list_all_valid_barcodes

        # specify information about R1 / R2 matches
        self._bcend = {
                'R1':self.bclen + len(self.downstream),
                'R2':self.bclen + len(self.upstream)
                }
        self._rcdownstream = dms_tools2.utils.reverseComplement(self.downstream)
        self._rcupstream = dms_tools2.utils.reverseComplement(self.upstream)
        self._matches = {'R1':{}, 'R2':{}} # saves match object by read length


    def parse(self, r1files, r2files=None):
        """Parses barcodes from files.

        Args:
            `r1file` (str or list)
                Name of R1 FASTQ file, or list of such files
                Can optionally be gzipped.
            `r2file` (`None`, str, or list)
                `None` of not using R2, otherwise like R1.

        Returns:
            The 2-tuple `(barcodes, fates)`. In this 2-tuple:

                - `barcodes` is a pandas DataFrame giving the
                  number of observations of each barcode. The
                  columns are named "barcode" and "count".

                - `fates` is a pandas DataFrame giving the
                  total number of reads with each fate. The
                  columns are named "fate" and "count".

                  - "valid barcode"

                  - "invalid barcode": not in our barcode whitelist

                  - "R1 / R2 disagree"

                  - "low quality barcode": sequencing quality low

                  - "unparseable barcode": invalid flanking sequences
                    or N in barcode.
        """
        if r2files is None:
            reads = ['R1']
        else:
            reads = ['R1', 'R2']

        if self.valid_barcodes and self.list_all_valid_barcodes:
            barcodes = {bc:0 for bc in self.valid_barcodes}
        else:
            barcodes = collections.defaultdict(int)

        fates = collections.defaultdict(int)

        for name, r1, r2, q1, q2, fail in \
                dms_tools2.utils.iteratePairedFASTQ(r1files, r2files):

            if fail and self.chastity_filter:
                fates['failed chastity filter'] += 1
                continue

            matches = {}
            for read, r in zip(reads, [r1, r2]):
                rlen = len(r)

                # get or build matcher for read of this length
                len_past_bc = rlen - self._bcend[read]
                if len_past_bc < 0:
                    raise ValueError(f"{read} too short: {rlen}")
                elif rlen in self._matches[read]:
                    matcher = self._matches[read][rlen]
                else:
                    if read == 'R1':
                        match_str = (
                                f'^({self._rcdownstream})'
                                f'{{s<={self.downstream_mismatch}}}' +
                                f'(?P<bc>N{{{self.bclen}}})' +
                                f'({self._rcupstream[ : len_past_bc]})' +
                                f'{{s<={self.upstream_mismatch}}}'
                                )
                    else:
                        assert read == 'R2'
                        match_str = (
                                f'^({self.upstream})' +
                                f'{{s<={self.upstream_mismatch}}}' +
                                f'(?P<bc>N{{{self.bclen}}})' +
                                f'({self.downstream[ : len_past_bc]})' +
                                f'{{s<={self.downstream_mismatch}}}'
                                )
                    matcher = regex.compile(
                            dms_tools2.pacbio.re_expandIUPAC(match_str),
                            flags=regex.BESTMATCH)
                    self._matches[read][rlen] = matcher

                m = matcher.match(r)
                if m:
                    matches[read] = m
                else:
                    break

            if len(matches) == len(reads):
                bc = {}
                bc_q = {}
                for read, q in zip(reads, [q1, q2]):
                    bc[read] = matches[read].group('bc')
                    bc_q[read] = numpy.array([
                                 ord(qi) - 33 for qi in 
                                 q[matches[read].start('bc') :
                                   matches[read].end('bc')]],
                                 dtype='int')
                if self.rc_barcode and 'R2' in reads:
                    bc['R2'] = dms_tools2.utils.reverseComplement(bc['R2'])
                    bc_q['R2'] = numpy.flip(bc_q['R2'], axis=0)
                elif 'R2' in reads:
                    bc['R1'] = dms_tools2.utils.reverseComplement(bc['R1'])
                    bc_q['R1'] = numpy.flip(bc_q['R1'], axis=0)
                if len(reads) == 1:
                    if (bc_q['R1'] >= self.minq).all():
                        if self.valid_barcodes and (
                                bc['R1'] not in self.valid_barcodes):
                            fates['invalid barcode'] += 1
                        else:
                            barcodes[bc['R1']] += 1
                            fates['valid barcode'] += 1
                    else:
                        fates['low quality barcode'] += 1
                else:
                    if bc['R1'] == bc['R2']:
                        if self.valid_barcodes and (
                                bc['R1'] not in self.valid_barcodes):
                            fates['invalid barcode'] += 1
                        elif (numpy.maximum(bc_q['R1'], bc_q['R2'])
                                >= self.minq).all():
                            barcodes[bc['R1']] += 1
                            fates['valid barcode'] += 1
                        else:
                            fates['low quality barcode'] += 1
                    else:
                        fates['R1 / R2 disagree'] += 1
            else:
                # invalid flanking sequence or N in barcode
                fates['unparseable barcode'] += 1

        barcodes = (pandas.DataFrame(
                        list(barcodes.items()),
                        columns=['barcode', 'count'])
                    .sort_values(['count', 'barcode'],
                                 ascending=[False, True])
                    .reset_index(drop=True)
                    )

        fates = (pandas.DataFrame(
                    list(fates.items()),
                    columns=['fate', 'count'])
                 .sort_values(['count', 'fate'],
                              ascending=[False, True])
                 .reset_index(drop=True)
                 )

        return (barcodes, fates)


class CodonVariantTable:
    """Associates barcodes with codon mutants of gene.

    Args:
        `barcode_variant_file` (str)
            CSV file giving barcodes and variants. Must have
            columns named "library", "barcode", "substitutions",
            (nucleotide mutations in 1, ... numbering in a format
            like "G301A A302T G856C"), and "variant_call_support"
            (sequences supporting barcode-variant call).
        `geneseq` (str)
            Sequence of protein-coding gene.

    Attributes:
        `sites` (list)
            List of all codon sites in 1, 2, ... numbering.
        `codons` (dict)
            `codons[r]` is wildtype codon at site `r`.
        `aas` (dict)
            `aas[r]` is wildtype amino acid at site `r`.
        `libraries` (list)
            List of all libraries in `barcode_variantfile`.
        `barcode_variant_df` (pandas DataFrame)
            Info about codon mutations parsed from `barcode_variantfile`.
        `variant_count_df` (pandas DataFrame or None)
            Initially `None`, but after data are added with
            :class:`CodonVariantTable.addSampleCounts`, holds
            counts of each variant for each sample. Differs from
            `barcode_variant_df` in that the former just holds
            initial barcode-variant data, whereas `variant_count_df`
            is updated with variant counts for samples.

    Initialize a :class:`CodonVariantTable`:

    >>> geneseq = 'ATGGGATGA'
    >>> variantfile = '_variantfile.csv'
    >>> with open(variantfile, 'w') as f:
    ...     _ = f.write(
    ...           'library,barcode,substitutions,variant_call_support\\n'
    ...           'lib_1,AAC,,2\\n'
    ...           'lib_1,GAT,G4C A6C,1\\n'
    ...           'lib_2,AAC,T2A G4T,2\\n'
    ...           'lib_2,CAT,A6C,3'
    ...           )
    >>> variants = CodonVariantTable(
    ...             barcode_variant_file=variantfile,
    ...             geneseq=geneseq
    ...             )
    >>> os.remove(variantfile)

    Check attributes of the :class:`CodonVariantTable`:

    >>> variants.sites
    [1, 2, 3]
    >>> variants.codons == {1:'ATG', 2:'GGA', 3:'TGA'}
    True
    >>> variants.aas == {1:'M', 2:'G', 3:'*'}
    True
    >>> variants.libraries
    ['lib_1', 'lib_2']
    >>> variants.valid_barcodes('lib_1') == {'AAC', 'GAT'}
    True
    >>> variants.valid_barcodes('lib_2') == {'AAC', 'CAT'}
    True
    >>> pandas.set_option('display.max_columns', 10)
    >>> pandas.set_option('display.width', 500)
    >>> variants.barcode_variant_df
      library barcode  variant_call_support codon_substitutions aa_substitutions  n_codon_substitutions  n_aa_substitutions
    0   lib_1     AAC                     2                                                           0                   0
    1   lib_1     GAT                     1             GGA2CGC              G2R                      1                   1
    2   lib_2     AAC                     2     ATG1AAG GGA2TGA          M1K G2*                      2                   2
    3   lib_2     CAT                     3             GGA2GGC                                       1                   0

    Initially we haven't added any barcode count information
    for any samples:

    >>> all(variants.samples(lib) == [] for lib in variants.libraries)
    True
    >>> variants.variant_count_df is None
    True

    Now we add barcode count information for a sample named
    "input" from library 1:

    >>> counts_lib1_input = pandas.DataFrame(
    ...         {'barcode':['AAC', 'GAT'],
    ...          'count'  :[  253,  1101]})
    >>> variants.addSampleCounts('lib_1', 'input', counts_lib1_input)
    >>> variants.variant_count_df
      barcode  count library sample  variant_call_support codon_substitutions aa_substitutions  n_codon_substitutions  n_aa_substitutions
    0     GAT   1101   lib_1  input                     1             GGA2CGC              G2R                      1                   1
    1     AAC    253   lib_1  input                     2                                                           0                   0

    We get an error if we try to add these same data again,
    as they are already added for that sample to that library:

    >>> variants.addSampleCounts('lib_1', 'input', counts_lib1_input)
    Traceback (most recent call last):
    ...
    ValueError: `library` lib_1 already has `sample` input

    But, we can add barcode counts for another
    sample (named "selected" in this case) to library 1:

    >>> counts_lib1_selected = pandas.DataFrame(
    ...         {'barcode':['AAC', 'GAT'],
    ...          'count'  :[  513,  401]})
    >>> variants.addSampleCounts('lib_1', 'selected', counts_lib1_selected)

    As well as barcode counts for the same two samples
    ("input" and "selected") to our other library (library 2):

    >>> counts_lib2_input = pandas.DataFrame(
    ...         {'barcode':['AAC', 'CAT'],
    ...          'count'  :[ 1253,  923]})
    >>> variants.addSampleCounts('lib_2', 'input', counts_lib2_input)
    >>> counts_lib2_selected = pandas.DataFrame(
    ...         {'barcode':['AAC', 'CAT'],
    ...          'count'  :[  113,  1200]})
    >>> variants.addSampleCounts('lib_2', 'selected', counts_lib2_selected)
    >>> variants.variant_count_df
      barcode  count library    sample  variant_call_support codon_substitutions aa_substitutions  n_codon_substitutions  n_aa_substitutions
    0     GAT   1101   lib_1     input                     1             GGA2CGC              G2R                      1                   1
    1     AAC    253   lib_1     input                     2                                                           0                   0
    2     AAC    513   lib_1  selected                     2                                                           0                   0
    3     GAT    401   lib_1  selected                     1             GGA2CGC              G2R                      1                   1
    4     AAC   1253   lib_2     input                     2     ATG1AAG GGA2TGA          M1K G2*                      2                   2
    5     CAT    923   lib_2     input                     3             GGA2GGC                                       1                   0
    6     CAT   1200   lib_2  selected                     3             GGA2GGC                                       1                   0
    7     AAC    113   lib_2  selected                     2     ATG1AAG GGA2TGA          M1K G2*                      2                   2
    """

    def __init__(self, *, barcode_variant_file, geneseq):
        """See main class doc string."""

        self.geneseq = geneseq.upper()
        if not re.match('^[ATGC]+$', self.geneseq):
            raise ValueError(f"invalid nucleotides in {self.geneseq}")
        if ((len(geneseq) % 3) != 0) or len(geneseq) == 0:
            raise ValueError(f"`geneseq` of invalid length {len(self.geneseq)}")
        self.sites = list(range(1, len(self.geneseq) // 3 + 1))
        self.codons = {r:self.geneseq[3 * (r - 1) : 3 * r] for r in self.sites}
        self.aas = {r:CODON_TO_AA[codon] for r, codon in self.codons.items()}

        df = pandas.read_csv(barcode_variant_file)
        required_cols = {'library', 'barcode',
                         'substitutions', 'variant_call_support'}
        if not set(df.columns).issuperset(required_cols):
            raise ValueError("`variantfile` does not have "
                             f"required columns {required_cols}")
        self.libraries = sorted(df.library.unique().tolist())
        self._valid_barcodes = {}
        for lib in self.libraries:
            barcodes = df.query('library == @lib').barcode
            if len(set(barcodes)) != len(barcodes):
                raise ValueError(f"duplicated barcodes for {lib}")
            self._valid_barcodes[lib] = set(barcodes)

        self._samples = {lib:[] for lib in self.libraries}
        self.variant_count_df = None

        self.barcode_variant_df = (
                df
                # info about codon and amino-acid substitutions
                .assign(codon_substitutions=
                            lambda x: x.substitutions
                                       .fillna('')
                                       .apply(self._ntToCodonMuts),
                        aa_substitutions=
                            lambda x: x.codon_substitutions
                                       .apply(self._codonToAAMuts),
                        n_codon_substitutions=
                            lambda x: x.codon_substitutions
                                       .str
                                       .split()
                                       .apply(len),
                        n_aa_substitutions=
                            lambda x: x.aa_substitutions
                                       .str
                                       .split()
                                       .apply(len)
                        )
                # we no longer need initial `substitutions` column
                .drop('substitutions', axis='columns')
                # sort to ensure consistent order
                .assign(library=lambda x:
                                pandas.Categorical(
                                    x.library,
                                    categories=self.libraries,
                                    ordered=True
                                    )
                        )
                .sort_values(['library', 'barcode'])
                .reset_index(drop=True)
                )


    def samples(self, library):
        """List of all samples for `library`.

        Args:
            `library` (str)
                Valid `library` for the :class:`CodonVariantTable`.

        Returns:
            List of all samples for which barcode counts have
            been added.
        """
        try:
            return self._samples[library]
        except KeyError:
            raise ValueError(f"invalid `library` {library}")


    def addSampleCounts(self, library, sample, barcodecounts):
        """Add variant counts for a sample to `variant_count_df`.

        Args:
            `library` (str)
                Valid `library` for the :class:`CodonVariantTable`.
            `sample` (str)
                Sample name, must **not** already be in
                :class:`CodonVariantTable.samples` for `library`.
            `barcodecounts` (pandas DataFrame)
                Gives counts for each variant by barcode. Must
                have columns named "barcode" and "count". The
                "barcode" column must contain all the barcodes
                in :class:`CodonVariantTable.valid_barcodes` for
                `library`. Such data frames are returned
                by :class:`IlluminaBarcodeParser.parse`.
        """
        if library not in self.libraries:
            raise ValueError(f"invalid library {library}")

        if sample in self.samples(library):
            raise ValueError(f"`library` {library} already "
                             f"has `sample` {sample}")

        req_cols = ['barcode', 'count']
        if not set(barcodecounts.columns).issuperset(set(req_cols)):
            raise ValueError(f"`barcodecounts` lacks columns {req_cols}")
        if len(barcodecounts) != len(set(barcodecounts.barcode.unique())):
            raise ValueError("`barcodecounts` has non-unique barcodes")
        if set(barcodecounts.barcode.unique()) != self.valid_barcodes(library):
            raise ValueError("barcodes in `barcodecounts` do not match "
                             f"those expected for `library` {library}")

        self._samples[library].append(sample)

        df = (barcodecounts
              [req_cols]
              .assign(library=library, sample=sample)
              .merge(self.barcode_variant_df,
                     how='inner',
                     on=['library', 'barcode'],
                     sort=False,
                     validate='one_to_one')
              )

        if self.variant_count_df is None:
            self.variant_count_df = df
        else:
            self.variant_count_df = pandas.concat(
                              [self.variant_count_df, df],
                              axis='index',
                              ignore_index=True,
                              sort=False
                              )

        # samples in order added after ordering by library, getting
        # unique ones as here: https://stackoverflow.com/a/39835527
        unique_samples = list(collections.OrderedDict.fromkeys(
                itertools.chain.from_iterable(
                    [self.samples(lib) for lib in self.libraries])
                ))

        # make library and sample categorical and sort
        self.variant_count_df = (
                self.variant_count_df
                .assign(library=lambda x:
                                pandas.Categorical(
                                    x.library,
                                    categories=self.libraries,
                                    ordered=True
                                    ),
                        sample=lambda x:
                               pandas.Categorical(
                                    x['sample'],
                                    categories=unique_samples,
                                    ordered=True
                                    )
                         )
                .sort_values(['library', 'sample', 'count'],
                             ascending=[True, True, False])
                .reset_index(drop=True)
                )


    def valid_barcodes(self, library):
        """Set of valid barcodes for `library`."""
        if library not in self.libraries:
            raise ValueError(f"invalid `library` {library}")
        else:
            return self._valid_barcodes[library]


    def plotNumCodonMutsByType(self, variant_type, *,
            libraries='all', samples='all', plotfile=None,
            orientation='h', widthscale=1, heightscale=1,
            min_support=1):
        """Nonsynonymous, synonymous, stop mutations per variant.

        Args:
            `variant_type` ("single" or "all")
                Include just single-codon mutants and wildtype,
                or include all mutants.
            Other args:
                Same meaning as for
                :class:`CodonVariantTable.plotnNumMutsHistogram`

        Returns:
            A `plotnine <https://plotnine.readthedocs.io/en/stable/>`_
        """
        df, nlibraries, nsamples = self._getPlotData(libraries,
                                                     samples,
                                                     min_support)

        if variant_type == 'single':
            df = df.query('n_codon_substitutions <= 1')
        elif variant_type != 'all':
            raise ValueError(f"invalid variant_type {variant_type}")

        if orientation == 'h':
            facet_str = 'sample ~ library'
            width = widthscale * (1 + 1.3 * nlibraries)
            height = heightscale * (1 + 1.2 * nsamples)
        elif orientation == 'v':
            facet_str = 'library ~ sample'
            width = widthscale * (1 + 1.2 * nsamples)
            height = heightscale * (1 + 1.2 * nlibraries)
        else:
            raise ValueError(f"invalid `orientation` {orientation}")

        if height > 3:
            ylabel = f'mutations per variant ({variant_type} mutants)'
        else:
            ylabel = f'mutations per variant\n({variant_type} mutants)'

        codon_mut_types = ['nonsynonymous', 'synonymous', 'stop']

        # mutations from stop to another amino-acid counted as nonsyn
        df = (df
              .assign(
                  synonymous=lambda x: x.n_codon_substitutions -
                                       x.n_aa_substitutions,
                stop=lambda x: x.aa_substitutions.str
                               .findall('[A-Z]\d+\*').apply(len),
                nonsynonymous=lambda x: x.n_codon_substitutions -
                                        x.synonymous - x.stop
                )
              .melt(id_vars=['library', 'sample'],
                    value_vars=codon_mut_types,
                    var_name='mutation_type',
                    value_name='number')
              .assign(
                  mutation_type=lambda x:
                                pandas.Categorical(
                                 x.mutation_type,
                                 categories=codon_mut_types,
                                 ordered=True),
                  nvariants=1
                  )
              .groupby(['library', 'sample', 'mutation_type'])
              .aggregate({'number':'sum', 'nvariants':'count'})
              .reset_index()
              .assign(number=lambda x: x.number / x.nvariants)
              )

        p = (ggplot(df, aes('mutation_type', 'number',
                            fill='mutation_type', label='number')) +
             geom_bar(stat='identity') +
             geom_text(size=8, va='bottom', format_string='{0:.3f}') +
             facet_grid(facet_str) +
             scale_y_continuous(name=ylabel,
                                expand=(0.03, 0, 0.12, 0)) +
             scale_fill_manual(COLOR_BLIND_PALETTE_GRAY[1 : ]) +
             theme(figure_size=(width, height),
                   axis_title_x=element_blank(),
                   axis_text_x=element_text(angle=90, size=11),
                   legend_position='none')
             )

        if plotfile:
            p.save(plotfile, height=height, width=width, verbose=False)

        return p


    def plotNumMutsHistogram(self, mut_type, *,
            libraries='all', samples='all', plotfile=None,
            orientation='h', widthscale=1, heightscale=1,
            min_support=1, max_muts=None):
        """Plot histograms of number of mutations per variant.

        Args:
            `mut_type` (str)
                Type of mutation to plot: "codon" or "aa".
            `libraries` (the str "all" or list)
                Which libraries do we plot? Set to "all" to plot
                all libraries plus sum over libraries, or a list
                of libraries to plot in indicated order.
            `samples` (the str "all", `None`, or list)
                Which samples do we plot? Set to "all" to plot all
                samples, `None` to plot counts of barcoded variants,
                or a list of samples to plot in indicated order.
            `plotfile` (`None` or str)
                Name of file to which we save plot.
            `orientation` (the str 'h' or 'v')
                Do we facet libraries horizontally or vertically?
            `widthscale` (float or int)
                Expand width of plot by this much.
            `heightscale` (float or int)
                Expand height of plot by this much.
            `min_support` (int)
                Only plot variants with at least this large
                of a variant call support.
            `max_muts` (int or `None`)
                In histogram, group together all variants with
                >= this many mutations; set to `None` for no
                cutoff.

        Returns:
            A `plotnine <https://plotnine.readthedocs.io/en/stable/>`_
            plot; can be displayed in a Jupyter notebook with `p.draw()`.
        """
        df, nlibraries, nsamples = self._getPlotData(libraries,
                                                     samples,
                                                     min_support)

        if mut_type == 'aa':
            mut_col = 'n_aa_substitutions'
            xlabel = 'amino-acid mutations'
        elif mut_type == 'codon':
            mut_col = 'n_codon_substitutions'
            xlabel = 'codon mutations'
        else:
            raise ValueError(f"invalid mut_type {mut_type}")

        df[mut_col] = numpy.clip(df[mut_col], None, max_muts)

        if orientation == 'h':
            facet_str = 'sample ~ library'
            width = widthscale * (1 + 1.2 * nlibraries)
            height = heightscale * (0.8 + 1.2 * nsamples)
        elif orientation == 'v':
            facet_str = 'library ~ sample'
            width = widthscale * (1 + 1.2 * nsamples)
            height = heightscale * (0.8 + 1.2 * nlibraries)
        else:
            raise ValueError(f"invalid `orientation` {orientation}")

        p = (ggplot(df, aes(mut_col)) +
             geom_histogram(binwidth=1, center=0) +
             facet_grid(facet_str) +
             xlab(xlabel) +
             theme(figure_size=(width, height))
             )

        if plotfile:
            p.save(plotfile, height=height, width=width, verbose=False)

        return p


    def _getPlotData(self, libraries, samples, min_support):
        """Gets data to plot from library and sample filters.

        Args:
            `libraries`, `samples`, `min_support` have meaning as
            for :class:`CodonVariantTable.plotNumMutsHistogram`.

        Returns:
            The 3-tuple `(df, nlibraries, nsamples)` where:

                - `df`: DataFrame with data to plot.

                - `nlibraries`: number of libraries being plotted.

                - `nsamples`: number of samples being plotted.
        """
        if samples is None:
            df = (self.barcode_variant_df
                  .assign(sample='barcoded variants')
                  )
        elif samples == 'all':
            if self.variant_count_df is None:
                raise ValueError('no samples have been added')
            df = self.variant_count_df
        elif isinstance(samples, list):
            all_samples = set(itertools.chain.from_iterable(
                    self.samples(lib) for lib in self.libraries))
            if not all_samples.issuperset(set(samples)):
                raise ValueError(f"invalid sample(s) in {samples}")
            if len(samples) != len(set(samples)):
                raise ValueError(f"duplicate samples in {samples}")
            df = (df
                  .query('sample in @samples')
                  .assign(sample=lambda x:
                                 pandas.Categorical(
                                   x['sample'],
                                   categories=samples,
                                   ordered=True)
                          )
                  )
        else:
            raise ValueError(f"invalid `samples` {samples}")

        df = df.query('variant_call_support >= @min_support')

        if not len(df):
            raise ValueError(f"no samples {samples}")
        else:
            nsamples = len(df['sample'].unique())

        if libraries == 'all':
            all_lib_name = 'all libraries'
            if all_lib_name in self.libraries:
                raise ValueError(f"library named {all_lib_name}")
            df = (pandas.concat([df, df.assign(library=all_lib_name)],
                                axis='index',
                                ignore_index=True,
                                sort=False)
                  .assign(library=lambda x:
                                  pandas.Categorical(
                                   x['library'],
                                   categories=self.libraries +
                                              [all_lib_name],
                                   ordered=True)
                          )
                  )
        elif isinstance(libraries, list):
            if not set(self.libraries).issuperset(set(libraries)):
                raise ValueError(f"invalid library in {libraries}")
            if len(libraries) != len(set(libraries)):
                raise ValueError(f"duplicate library in {libraries}")
            df = (df
                  .query('library in @libraries')
                  .assign(library=lambda x:
                                  pandas.Categorical(
                                   x['library'],
                                   categories=libraries,
                                   ordered=True)
                          )
                  )
        else:
            raise ValueError(f"invalid `libraries` {libraries}")
        if not len(df):
            raise ValueError(f"no libraries {libraries}")
        else:
            nlibraries = len(df['library'].unique())

        return (df, nlibraries, nsamples)


    def _codonToAAMuts(self, codon_mut_str):
        """Converts string of codon mutations to amino-acid mutations.

        Args:
            `codon_mut_str` (str)
                Codon mutations, delimited by a space and in
                1, 2, ... numbering.

        Returns:
            String with amino acid mutations in 1, 2, ... numbering.

        >>> geneseq = 'ATGGGATGA'
        >>> with tempfile.NamedTemporaryFile(mode='w') as f:
        ...     _ = f.write('library,barcode,substitutions,variant_call_support')
        ...     f.flush()
        ...     variants = CodonVariantTable(
        ...                 barcode_variant_file=f.name,
        ...                 geneseq=geneseq
        ...                 )
        >>> variants._codonToAAMuts('ATG1GTG GGA2GGC TGA3AGA')
        'M1V *3R'
        """
        aa_muts = {}
        for mut in codon_mut_str.upper().split():
            m = re.match('^(?P<wt>[ATGC]{3})(?P<r>\d+)(?P<mut>[ATGC]{3})$',
                         mut)
            if not m:
                raise ValueError(f"invalid codon mutation {mut}")
            r = int(m.group('r'))
            if r in aa_muts:
                raise ValueError(f"duplicate codon mutation for {r}")
            wt_codon = m.group('wt')
            if self.geneseq[3 * (r - 1) : 3 * r] != wt_codon:
                raise ValueError(f"invalid wildtype codon in {mut}")
            mut_codon = m.group('mut')
            if wt_codon == mut_codon:
                raise ValueError(f"invalid mutation {mut}")
            wt_aa = dms_tools2.CODON_TO_AA[wt_codon]
            assert wt_aa == self.aas[r]
            mut_aa = dms_tools2.CODON_TO_AA[mut_codon]
            if wt_aa != mut_aa:
                aa_muts[r] = f"{wt_aa}{r}{mut_aa}"

        return ' '.join([mut_str for r, mut_str in sorted(aa_muts.items())])


    def _ntToCodonMuts(self, nt_mut_str):
        """Converts string of nucleotide mutations to codon mutations.

        Args:
            `nt_mut_str` (str)
                Nucleotide mutations, delimited by a space and in
                1, 2, ... numbering.

        Returns:
            String with codon mutations in 1, 2, ... numbering of
            codon sites.

        >>> geneseq = 'ATGGGATGA'
        >>> with tempfile.NamedTemporaryFile(mode='w') as f:
        ...     _ = f.write('library,barcode,substitutions,variant_call_support')
        ...     f.flush()
        ...     variants = CodonVariantTable(
        ...                 barcode_variant_file=f.name,
        ...                 geneseq=geneseq
        ...                 )
        >>> variants._ntToCodonMuts('A1G G4C A6T')
        'ATG1GTG GGA2CGT'
        >>> variants._ntToCodonMuts('A1G G4C G6T')
        Traceback (most recent call last):
        ...
        ValueError: nucleotide 6 should be A not G
        """
        mut_codons = collections.defaultdict(set)
        for mut in nt_mut_str.upper().split():
            m = re.match('^(?P<wt>[ATCG])(?P<i>\d+)(?P<mut>[ATCG])$', mut)
            if not m:
                raise ValueError(f"invalid mutation {mut}")
            wt_nt = m.group('wt')
            i = int(m.group('i'))
            mut_nt = m.group('mut')
            if wt_nt == mut_nt:
                raise ValueError(f"invalid mutation {mut}")
            if i > len(self.geneseq):
                raise ValueError(f"invalid nucleotide site {i}")
            if self.geneseq[i - 1] != wt_nt:
                raise ValueError(f"nucleotide {i} should be "
                                 f"{self.geneseq[i - 1]} not {wt_nt}")
            icodon = (i - 1) // 3 + 1
            i_nt = (i - 1) % 3
            assert self.codons[icodon][i_nt] == wt_nt
            if i_nt in mut_codons[icodon]:
                raise ValueError(f"duplicate mutations {i_nt} in {icodon}")
            mut_codons[icodon].add((i_nt, mut_nt))

        codon_mut_list = []
        for r, r_muts in sorted(mut_codons.items()):
            wt_codon = self.codons[r]
            mut_codon = list(wt_codon)
            for i, mut_nt in r_muts:
                mut_codon[i] = mut_nt
            codon_mut_list.append(f"{wt_codon}{r}{''.join(mut_codon)}")

        return ' '.join(codon_mut_list)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
