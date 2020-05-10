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
import Bio.SeqUtils.ProtParamData

# use plotnine for plotting
from plotnine import *

import dms_tools2.plot
from dms_tools2.plot import COLOR_BLIND_PALETTE_GRAY
import dms_tools2.utils
import dms_tools2.pacbio
from dms_tools2 import CODON_TO_AA, CODONS, AAS_WITHSTOP, AA_TO_CODONS

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
        if isinstance(barcodes, collections.abc.Iterable):
            barcodes = pandas.Series(barcodes)
        else:
            raise TypeError(f"`barcodes` invalid type {type(barcodes)}")
    barcodes = barcodes.str.encode('utf-8')

    counts = collections.Counter(barcodes)

    groups = []
    for group in _umi_clusterer(counts, threshold):
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
    >>> pandas.set_option('display.max_columns', 10)
    >>> pandas.set_option('display.width', 500)
    >>> consensus
      library barcode substitutions  insertions deletions  variant_call_support
    0      s1      AG         [A2C]          []        []                     2
    1      s1      TA         [G3A]  [ins4len3]        []                     1
    2      s2      GG            []          []        []                     2
    3      s2      TA         [T6C]          []        []                     3
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
                    m = re.search(r'(\-{0,1}\d+)', mut)
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

        5'-[R2_start]-upstream-barcode-downstream-[R1_start]-3'

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
        `bclen` (int or `None`)
            Length of the barcode, or `None` if length is to be
            determined from `valid_barcodes` argument.
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
    ...         '????????????\\n')
    ...     _ = f2.write(
    ...         '@valid_CGTA_barcode_partial_flanking_region\\n'
    ...         'ACATGATACG\\n'
    ...         '+\\n'
    ...         '??????????\\n')
    ...
    ...     # valid GCCG barcode, extended flanking regions
    ...     _ = f1.write(
    ...         '@valid_GCCG_barcode_extended_flanking_region\\n'
    ...         'AGTCGCCGTCATGTTAC\\n'
    ...         '+\\n'
    ...         '?????????????????\\n')
    ...     _ = f2.write(
    ...         '@valid_GCCG_barcode_extended_flanking_region\\n'
    ...         'ACATGACGGCGACTGAC\\n'
    ...         '+\\n'
    ...         '?????????????????\\n')
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

    def __init__(self, *, bclen=None,
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
        if self.valid_barcodes is not None:
            self.valid_barcodes = set(self.valid_barcodes)
            if len(self.valid_barcodes) < 1:
                raise ValueError('empty list for `valid_barcodes`')
            if self.bclen is None:
                self.bclen = len(list(self.valid_barcodes)[0])
            if any(len(bc) != self.bclen for bc in self.valid_barcodes):
                raise ValueError('`valid_barcodes` not all valid length')
        elif self.bclen is None:
            raise ValueError('must specify `bclen` or `valid_barcodes`')
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
            `r1files` (str or list)
                Name of R1 FASTQ file, or list of such files
                Can optionally be gzipped.
            `r2files` (`None`, str, or list)
                `None` or empty list if not using R2, otherwise like R1.

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
        if not r2files:
            reads = ['R1']
            r2files = None
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


def tidy_split(df, column, sep=' ', keep=False):
    """
    Split values of a column and expand so new DataFrame has one split
    value per row. Filters rows where the column is missing.

    Taken from https://stackoverflow.com/a/39946744

    Args:
        df : pandas DataFrame
            dataframe with the column to split and expand
        column : str
            the column to split and expand
        sep : str
            the string used to split the column's values
        keep : bool
            whether to retain the presplit value as it's own row

    Returns:
        pandas DataFrame
            Returns a dataframe with the same columns as `df`.
    """
    indexes = list()
    new_values = list()
    df = df.dropna(subset=[column])
    for i, presplit in enumerate(df[column].astype(str)):
        values = presplit.split(sep)
        if keep and len(values) > 1:
            indexes.append(i)
            new_values.append(presplit)
        for value in values:
            indexes.append(i)
            new_values.append(value)
    new_df = df.iloc[indexes, :].copy()
    new_df[column] = new_values
    return new_df


if __name__ == '__main__':
    import doctest
    doctest.testmod()
