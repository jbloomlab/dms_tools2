"""
===================
pacbio
===================

Tools for processing PacBio sequencing data.
"""


import os
import re
import io
import math
import subprocess
import collections
import tempfile

import numpy
import pandas
import pysam

# import dms_tools2.plot to set plotting contexts / themes
import dms_tools2.plot
from dms_tools2.plot import COLOR_BLIND_PALETTE
from dms_tools2.plot import COLOR_BLIND_PALETTE_GRAY
import matplotlib.pyplot as plt
import seaborn
from plotnine import *


class CCS:
    """Class to handle results of ``ccs``.

    Holds results of PacBio ``ccs``.
    Has been tested on output of ``ccs`` version 3.0.0.

    This class reads all data into memory, and so you
    may need a lot of RAM if `bamfile` is large.

    Args:
        `samplename` (str)
            Sample or sequencing run
        `bamfile` (str)
            BAM file created by ``ccs``
        `reportfile` (str or `None`)
            Report file created by ``ccs``, or
            `None` if you have no reports.

    Attributes:
        `samplename` (str)
            Name set at initialization
        `bamfile` (str)
            ``ccs`` BAM file set at initialization
        `reportfile` (str or `None`)
            ``ccs`` report file set at initialization
        `zmw_report` (pandas.DataFrame or `None`):
            ZMW stats in `reportfile`, or `None` if no
            `reportfile`. Columns are *status*, *number*,
            *percent*, and *fraction*.
        `subread_report` (pandas.DataFrame or `None`)
            Like `zmw_report` but for subreads.
        `df` (pandas.DataFrame)
            The CCSs in `bamfile`. Each row is a different CCS
            On creation, there will be the following columns (you
            can modify to add more): 

              - "name": the name of the CCS
              - "samplename": the sample as set via `samplename`
              - "CCS": the circular consensus sequence
              - "CCS_qvals": the Q-values as a numpy array
              - "passes": the number of passes of the CCS
              - "CCS_accuracy": the accuracy of the CCS
              - "CCS_length": the length of the CCS

    Here is an example.

    First, define the sequences, quality scores,
    and names for 3 example sequences. The names indicate
    the barcodes, the accuracy of the barcode, and the polarity.
    Two of the sequences have the desired termini and
    a barcode. The other does not. Note that the second
    sequence has an extra nucleotide at each end, this
    will turn out to be fine with the `match_str` we write.
    The second sequence is also reverse complemented:

    >>> termini5 = 'ACG'
    >>> termini3 = 'CTT'
    >>> ccs_seqs = [
    ...         {'name':'barcoded_TTC_0.999_plus',
    ...          'seq':termini5 + 'TTC' + 'ACG' + termini3,
    ...          'qvals':'?' * 12,
    ...         },
    ...         {'name':'barcoded_AGA_0.995_minus',
    ...          'seq':dms_tools2.utils.reverseComplement(
    ...                'T' + termini5 + 'AGA' + 'GCA' + termini3 + 'A'),
    ...          'qvals':''.join(reversed('?' * 4 + '5?9' + '?' * 7)),
    ...         },
    ...         {'name':'invalid',
    ...          'seq':'GGG' + 'CAT' + 'GCA' + termini3,
    ...          'qvals':'?' * 12,
    ...         }
    ...         ]
    >>> for iccs in ccs_seqs:
    ...     iccs['accuracy'] = qvalsToAccuracy(iccs['qvals'], encoding='sanger')

    Now place these in a block of text that meets the
    `CCS SAM specification <https://github.com/PacificBiosciences/unanimity/blob/develop/doc/PBCCS.md>`_:

    >>> sam_template = '\\t'.join([
    ...        '{0[name]}',
    ...        '4', '*', '0', '255', '*', '*', '0', '0',
    ...        '{0[seq]}',
    ...        '{0[qvals]}',
    ...        'np:i:6',
    ...        'rq:f:{0[accuracy]}',
    ...        ])
    >>> samtext = '\\n'.join([sam_template.format(iccs) for
    ...                      iccs in ccs_seqs])

    Create small SAM file with these sequences, then
    convert to BAM file used to initialize a `CCS` object
    (note this requires ``samtools`` to be installed):

    >>> samfile = '_temp.sam'
    >>> bamfile = '_temp.bam'
    >>> with open(samfile, 'w') as f:
    ...     _ = f.write(samtext)
    >>> _ = subprocess.check_call(['samtools', 'view',
    ...         '-b', '-o', bamfile, samfile])
    >>> ccs = CCS('test', bamfile, None)
    >>> os.remove(samfile)
    >>> os.remove(bamfile)

    Check `ccs.df` has correct names, samplename, CCS sequences,
    and columns:

    >>> set(ccs.df.name) == {s['name'] for s in ccs_seqs}
    True
    >>> all(ccs.df.samplename == 'test')
    True
    >>> set(ccs.df.CCS) == {s['seq'] for s in ccs_seqs}
    True
    >>> set(ccs.df.columns) == {'CCS', 'CCS_qvals', 'name',
    ...         'passes', 'CCS_accuracy', 'CCS_length', 'samplename'}
    True

    Match sequences that have expected termini
    and define barcode and read in these:

    >>> match_str = (termini5 + '(?P<barcode>N{3})' +
    ...         '(?P<read>N+)' + termini3)
    >>> ccs.df = matchSeqs(ccs.df, match_str, 'CCS', 'barcoded')

    This matching add new columns to the new `ccs.df`:

    >>> set(ccs.df.columns) >= {'barcode', 'barcode_qvals',
    ...         'barcode_accuracy', 'read', 'read_qvals',
    ...         'read_accuracy', 'barcoded', 'barcoded_polarity'}
    True

    Now make sure `df` indicates that the correct sequences
    are barcoded, and that they have the correct barcodes:

    >>> bc_names = sorted([s['name'] for s in ccs_seqs if
    ...         'barcoded' in s['name']])
    >>> ccs.df = ccs.df.sort_values('barcode')
    >>> (ccs.df.query('barcoded').name == bc_names).all()
    True
    >>> barcodes = [x.split('_')[1] for x in bc_names]
    >>> (ccs.df.query('barcoded').barcode == barcodes).all()
    True
    >>> (ccs.df.query('not barcoded').barcode == ['']).all()
    True
    >>> barcode_accuracies = [float(x.split('_')[2]) for x in bc_names]
    >>> numpy.allclose(ccs.df.query('barcoded').barcode_accuracy,
    ...     barcode_accuracies, atol=1e-4)
    True
    >>> numpy.allclose(ccs.df.query('barcoded').barcode_accuracy,
    ...         [qvalsToAccuracy(qvals) for qvals in
    ...         ccs.df.query('barcoded').barcode_qvals])
    True
    >>> numpy.allclose(ccs.df.query('not barcoded').barcode_accuracy,
    ...     -1, atol=1e-4)
    True
    >>> barcoded_polarity = [{'plus':1, 'minus':-1}[x.split('_')[3]]
    ...         for x in bc_names]
    >>> (ccs.df.query('barcoded').barcoded_polarity == barcoded_polarity).all()
    True

    """

    def __init__(self, samplename, bamfile, reportfile):
        """See main class doc string."""
        self.samplename = samplename

        assert os.path.isfile(bamfile), "can't find {0}".format(bamfile)
        self.bamfile = bamfile

        self.reportfile = reportfile
        if self.reportfile is None:
            self.zmw_report = None
            self.subread_report = None
        else:
            assert os.path.isfile(reportfile), \
                    "can't find {0}".format(reportfile)
            # set `zmw_report` and `subread_report`
            self._parse_report()

        self._build_df_from_bamfile()


    def plotColCorrs(self, plotfile, lower_filter=None,
            cols=['passes', 'CCS_accuracy', 'CCS_length'],
            title=True):
        """Plots correlation among CCS columns in `df`.

        The plot shows the distribution of each variable
        as well as all pairwise correlations.

        Args:
            `plotfile` (str)
                Name of created plot.
            `lower_filter` (`None` or str)
                Can specify boolean column in `df`. In this
                case, on lower diagonal only plot data
                for which this column is `True`.
            `cols` (list)
                List of variables to plot. There must be a
                column for each in `df`.
            `title` (bool or str)
                If `False`, no title. If `True`, make
                `samplename` the title. If a string, make
                that the title.
        """
        assert set(cols) <= set(self.df.columns)

        if lower_filter is not None:
            assert lower_filter in self.df.columns, \
                    "No `lower_filter` column {0}".format(lower_filter)
            assert self.df[lower_filter].dtype == 'bool', \
                    "`lower_filter` not boolean column"
            filter_indices = self.df.query(lower_filter).index
            color_all = COLOR_BLIND_PALETTE_GRAY[0]
            color_filter = COLOR_BLIND_PALETTE_GRAY[1]
        else:
            color_all = COLOR_BLIND_PALETTE[0]

        def hist1d(x, color, **kwargs):
            """1D histogram for diagonal elements."""
            bins=dms_tools2.plot.hist_bins_intsafe(x,
                    shrink_threshold=50)
            plt.hist(x, color=color_all, bins=bins, **kwargs)
            if lower_filter:
                plt.hist(x.ix[filter_indices], color=color_filter,
                         bins=bins, **kwargs)

        def hist2d(x, y, color, filterdata, **kwargs):
            """2D histogram for off-diagonal elements."""
            bins = [dms_tools2.plot.hist_bins_intsafe(a,
                    shrink_threshold=50) for a in [x, y]]
            if filterdata:
                color = color_filter
                x = x.ix[filter_indices]
                y = y.ix[filter_indices]
            else:
                color = color_all
            cmap = dms_tools2.plot.from_white_cmap(color)
            plt.hist2d(x, y, bins=bins, cmap=cmap, **kwargs)

        g = (dms_tools2.plot.AugmentedPairGrid(self.df, vars=cols,
                diag_sharey=False, size=3)
             .map_diag(hist1d)
             .map_upper(hist2d, filterdata=False)
             .map_lower(hist2d, filterdata=(lower_filter is not None))
             .ax_lims_clip_outliers()
             )

        if lower_filter is not None:
            label_order = ['all CCSs\n({0})'.format(
                                dms_tools2.plot.latexSciNot(len(self.df))),
                           '{0} CCSs\n({1})'.format(lower_filter,
                                dms_tools2.plot.latexSciNot(
                                        len(self.df.query(lower_filter)))),
                          ]
            label_data = {lab:plt.Line2D([0], [0], color=c, lw=10,
                                         solid_capstyle='butt')
                          for (lab, c) in zip(label_order,
                                [color_all, color_filter])
                          }
            g.add_legend(label_data, label_order=label_order,
                         labelspacing=2, handlelength=1.5)

        if title:
            if not isinstance(title, str):
                title = self.samplename
            g.fig.suptitle(title, va='bottom')

        g.savefig(plotfile)
        plt.close()

    def _parse_report(self):
        """Set `zmw_report` and `subread_report` using `reportfile`."""
        # match reports made by ccs 3.0.0
        reportmatch = re.compile('^ZMW Yield\n(?P<zmw>(.+\n)+)\n\n'
                            'Subread Yield\n(?P<subread>(.+\n)+)$')

        with open(self.reportfile) as f:
            report = f.read()
        m = reportmatch.search(report)
        assert m, "Cannot match {0}\n\n{1}".format(
                self.reportfile, report)

        for read_type in ['zmw', 'subread']:
            df = (pandas.read_csv(
                        io.StringIO(m.group(read_type)),
                        names=['status', 'number', 'percent']
                        )
                  .assign(fraction=lambda x: 
                        x.percent.str.slice(None, -1)
                        .astype('float') / 100)
                  )
            setattr(self, read_type + '_report', df)


    def _build_df_from_bamfile(self):
        """Builds `df` from `bamfile`."""
        # read into dictionary
        d = collections.defaultdict(list)
        for s in pysam.AlignmentFile(self.bamfile, 'rb',
                check_sq=False):
            d['CCS'].append(s.query_sequence)
            d['CCS_qvals'].append(numpy.asarray(s.query_qualities,
                                                dtype='int'))
            d['name'].append(s.query_name)
            d['passes'].append(s.get_tag('np'))
            d['CCS_accuracy'].append(s.get_tag('rq'))
            d['CCS_length'].append(s.query_length)
            d['samplename'].append(self.samplename)

        # create data frame
        self.df = pandas.DataFrame(d)

        # some checks on `df`
        assert self.df.name.size == self.df.name.unique().size,\
                "non-unique names for {0}".format(self.name)
        assert (self.df.CCS_length == self.df.CCS.apply(len)).all(),\
                "CCS not correct length"
        assert (self.df.CCS_length == self.df.CCS_qvals.apply(len)).all(),\
                "qvals not correct length"


def matchSeqs(df, match_str, col_to_match, match_col, *,
        add_polarity=True, add_group_cols=True,
        add_accuracy=True, add_qvals=True,
        expandIUPAC=True, overwrite=False):
    """Identify sequences in a dataframe that match a specific pattern.

    Args:
        `df` (pandas DataFrame)
            Data frame with column holding sequences to match.
        `match_str` (str)
            A string that can be passed to `re.compile` that gives
            the pattern that we are looking for, with target 
            subsequences as named groups. See also the `expandIUPAC`
            parameter, which simplifies writing `match_str`.
        `col_to_match` (str)
            Name of column in `df` that contains the sequences
            to match.
        `match_col` (str)
            Name of column added to `df`. Elements of columns are
            `True` if `col_to_match` matches `match_str` for that
            row, and `False` otherwise.
        `add_polarity` (bool)
            Add a column specifying the polarity of the match?
        `add_group_cols` (bool)
            Add columns with the sequence of every group in
            `match_str`?
        `add_accuracy` (bool)
            For each group in the match, add a column giving
            the accuracy of that group's sequence? Only used
            if `add_group_cols` is `True`.
        `add_qvals` (bool)
            For each group in the match, add a column giving
            the Q values for that group's sequence? Only used if
            `add_group_cols` is `True`.
        `expandIUPAC` (bool)
            Use `IUPAC code <https://en.wikipedia.org/wiki/Nucleic_acid_notation>`_
            to expand ambiguous nucleotides (e.g., "N") by passing
            `match_str` through the :meth:`re_expandIUPAC` function.
        `overwrite` (bool)
            If `True`, we overwrite any existing columns to
            be created that already exist. If `False`, raise
            an error if any of the columns already exist.

    Returns:
        A **copy** of `df` with new columns added. The exact columns
        to add are specified by the calling arguments. Specifically:

            - We always add a column with the name given by `match_col`
              that is `True` if there was a match and `False` otherwise.

            - If `add_polarity` is `True`, add a column that is
              `match_col` suffixed by "_polarity" which is 1 if
              the match is directly to the sequence in `col_to_match`,
              and -1 if it is to the reverse complement of this sequence.
              The value is 0 if there is no match.

            - If `add_group_cols` is `True`, then for each group
              in `match_str` specified using the `re` group naming
              syntax, add a column with that group name that
              gives the sequence matching that group. These
              sequences are empty strings if there is no match.
              These added sequences are in the polarity of the
              match, so if the sequence in `match_col` has
              to be reverse complemented for a match, then these
              sequences will be the reverse complement that matches.
              Additionally, when `add_group_cols` is True:

                - If `add_accuracy` is `True`, we also add a column
                  suffixed by "_accuracy" that gives the
                  accuracy of that group as computed from the Q-values.
                  The value -1 if there is match for that row. Adding
                  accuracy requires a colum in `df` with the name
                  given by `match_col` suffixed by "_qvals."

                - If `add_qvals` is `True`, we also add a column 
                  suffixed by "_qvals" that gives the Q-values
                  for that sequence. Adding these Q-values requires
                  that there by a column in `df` with the name given by
                  `match_col` suffixed by "_qvals". The Q-values are
                  in the form of a numpy array, or an empty numpy array
                  if there is no match for that row.
              
    See the docs for :class:`CCS` for example use of this function."""

    assert col_to_match in df.columns, \
            "`df` lacks `col_to_match` column {0}".format(col_to_match)

    if expandIUPAC:
        match_str = re_expandIUPAC(match_str)
    matcher = re.compile(match_str)

    newcols = [match_col]
    if add_polarity:
        polarity_col = match_col + '_polarity'
        newcols.append(polarity_col)

    if add_group_cols:
        groupnames = list(matcher.groupindex.keys())
        if len(set(groupnames)) != len(groupnames):
            raise ValueError("duplicate group names in {0}"
                             .format(match_str))
        newcols += groupnames
        if add_accuracy:
            newcols += [g + '_accuracy' for g in groupnames]
        if add_qvals:
            newcols += [g + '_qvals' for g in groupnames]
        if add_accuracy or add_qvals:
            match_qvals_col = col_to_match + '_qvals'
            if match_qvals_col not in df.columns:
                raise ValueError("To use `add_accuracy` or "
                        "`add_qvals`, you need a column in `df` "
                        "named {0}".format(match_qvals_col))
    else:
        groupnames = []

    # make sure created columns don't already exist
    dup_cols = set(newcols).intersection(set(df.columns))
    if not overwrite and dup_cols:
        raise ValueError("`df` already contains some of the "
                "columns that we are supposed to add:\n{0}"
                .format(dup_cols))

    # look for matches for each row
    match_d = {c:[] for c in newcols}
    for tup in df.itertuples():
        s = getattr(tup, col_to_match)
        m = matcher.search(s)
        if add_accuracy or add_qvals:
            qs = getattr(tup, match_qvals_col)
        if m:
            polarity = 1
        else:
            m = matcher.search(dms_tools2.utils.reverseComplement(s))
            polarity = -1
            if add_accuracy or add_qvals:
                qs = numpy.flip(qs, axis=0)
        if m:
            match_d[match_col].append(True)
            if add_polarity:
                match_d[polarity_col].append(polarity)
            for g in groupnames:
                match_d[g].append(m.group(g))
                if add_qvals:
                    match_d[g + '_qvals'].append(qs[m.start(g) : m.end(g)])
                if add_accuracy:
                    match_d[g + '_accuracy'].append(qvalsToAccuracy(
                            qs[m.start(g) : m.end(g)]))
        else:
            match_d[match_col].append(False)
            if add_polarity:
                match_d[polarity_col].append(0)
            for g in groupnames:
                match_d[g].append('')
                if add_qvals:
                    match_d[g + '_qvals'].append(numpy.array([], dtype='int'))
                if add_accuracy:
                    match_d[g + '_accuracy'].append(-1)

    # set index to make sure matches `df`
    indexname = df.index.name
    assert indexname not in match_d
    match_d[indexname] = df.index.tolist()
    if (not overwrite) and dup_cols:
        raise ValueError("overwriting columns")
    return pandas.concat(
            [df.drop(dup_cols, axis=1),
                pandas.DataFrame(match_d).set_index(indexname),
            ],
            axis=1)


def alignSeqs(df, mapper, query_col, aligned_col, *,
        add_alignment=True, add_target=True, add_cigar=True,
        add_endtoend=True, add_n_additional=True, overwrite=True,
        paf_file=None):
    """Align sequences in a dataframe to target sequence(s).

    Arguments:
        `df` (pandas DataFrame)
            Data frame in which one column holds sequences to match.
        `mapper` (:py:mod:`dms_tools2.minimap2.Mapper`)
            Align using the :py:mod:`dms_tools2.minimap2.Mapper.map`
            function of `mapper`. Target sequence(s) to which
            we align are specified when initializing `mapper`.
        `query_col` (str)
            Name of column in `df` with query sequences to align.
        `aligned_col` (str)
            Name of column added to `df`. Elements of column are
            `True` if `query_col` aligns, and `False` otherwise.
        `add_alignment` (bool)
            Add column with the :py:mod:`dms_tools2.minimap2.Alignment`.
        `add_target` (bool)
            Add column giving target (reference) to which sequence
            aligns.
        `add_cigar` (bool)
            Add column with the CIGAR string in the long format
            `described here <https://github.com/lh3/minimap2#cs>`_.
        `add_endtoend` (bool)
            Add column specifying whether the alignment is
            end-to-end in both query and target.
        `add_n_additional` (bool)
            Add column specifying the number of additional
            alignments.
        `overwrite` (bool)
            If `True`, we overwrite any existing columns to
            be created that already exist. If `False`, raise
            an error if any of the columns already exist.
        `paf_file` (`None` or str)
            If a str, is the name of the PAF file created
            by `mapper` (see `outfile` argument of
            :py:mod:`dms_tools2.minimap2.Mapper.map`) Otherwise
            this file is not saved.

    Returns:
        A **copy** of `df` with new columns added. The exact
        columns to add are specified by the calling arguments.
        Specifically:

            - We always add a column with the name given by
              `aligned_col` that is `True` if there was an
              alignment and `False` otherwise.
              
            - If `add_alignment` is `True`, add column named
              `aligned_col` suffixed by "_alignment" that gives
              the alignment as a :py:mod:`dms_tools2.minimap2.Alignment`
              object, or `None` if there is no alignment. Note that
              if there are multiple alignments, then this is the
              "best" alignment, and the remaining alignments are in
              the :py:mod:`dms_tools2.minimap2.Alignment.additional`
              attribute.

            - If `add_target` is `True`, add column named
              `aligned_col` suffixed by "_target" that gives
              the target to which the sequence aligns in the
              "best" alignment, or an empty string if no alignment.

            - If `add_cigar` is `True`, add column named
              `aligned_col` suffixed by "_cigar" with the CIGAR
              string (`long format <https://github.com/lh3/minimap2#cs>`_)
              for the "best" alignment, or an empty string if there
              is no alignment.

            - If `add_endtoend` is `True`, add column named
              `aligned_col` suffixed by "_endtoend" that is `True`
              if the "best" alignment goes from end-to-end in
              both the query and target. There can still be gaps in
              an end-to-end alignment, but they are internal.

            - If `add_n_additional` is `True`, add column
              named `aligned_col` suffixed by "_n_additional" that
              gives the number of additional alignments (in
              :py:mod:`dms_tools2.minimap2.Alignment.additional`),
              or -1 if there are no additional alignments.
    """
    assert query_col in df.columns, "no `query_col` {0}".format(query_col)

    newcols = [aligned_col]
    if add_alignment:
        alignment_col = aligned_col + '_alignment'
        newcols.append(alignment_col)
    if add_target:
        target_col = aligned_col + '_target'
        newcols.append(target_col)
    if add_cigar:
        cigar_col = aligned_col + '_cigar'
        newcols.append(cigar_col)
    if add_endtoend:
        endtoend_col = aligned_col + '_endtoend'
        newcols.append(endtoend_col)
    if add_n_additional:
        n_additional_col = aligned_col + '_n_additional'
        newcols.append(n_additional_col)

    dup_cols = set(newcols).intersection(set(df.columns))
    if (not overwrite) and dup_cols:
        raise ValueError("`df` already contains these columns:\n{0}"
                         .format(dup_cols))

    assert len(df.name) == len(df.name.unique()), \
            "`name` in `df` not unique"
    with tempfile.NamedTemporaryFile(mode='w') as queryfile:
        queryfile.write('\n'.join([
                        '>{0}\n{1}'.format(*tup) for tup in
                        df.query('{0} != ""'.format(query_col))
                            [['name', query_col]]
                            .itertuples(index=False, name=False)
                        ]))
        map_dict = mapper.map(queryfile.name, outfile=paf_file)

    align_d = {c:[] for c in newcols}
    for name in df.name:
        if name in map_dict:
            a = map_dict[name]
            assert a.strand == 1, "method does not handle - polarity"
            align_d[aligned_col].append(True)
            if add_alignment:
                align_d[alignment_col].append(a)
            if add_target:
                align_d[target_col].append(a.target)
            if add_cigar:
                align_d[cigar_col].append(a.cigar_str)
            if add_endtoend:
                if ((a.r_st == 0) and (a.q_st == 0) and
                        (a.q_en == a.q_len) and
                        (a.r_en == len(mapper.targetseqs[a.target]))):
                    align_d[endtoend_col].append(True)
                else:
                    align_d[endtoend_col].append(False)
            if add_n_additional:
                align_d[n_additional_col].append(len(a.additional))
        else:
            align_d[aligned_col].append(False)
            if add_alignment:
                align_d[alignment_col].append(None)
            if add_target:
                align_d[target_col].append('')
            if add_cigar:
                align_d[cigar_col].append('')
            if add_endtoend:
                align_d[endtoend_col].append(False)
            if add_n_additional:
                align_d[n_additional_col].append(-1)

    # set index to make sure matches `df`
    index_name = df.index.name
    assert index_name not in align_d
    align_d[index_name] = df.index.tolist()
    if (not overwrite) and dup_cols:
        raise ValueError("overwriting columns")
    return pandas.concat(
            [df.drop(dup_cols, axis=1),
                pandas.DataFrame(align_d).set_index(index_name),
            ],
            axis=1)


def qvalsToAccuracy(qvals, encoding='numbers'):
    """Converts set of quality scores into average accuracy.

    Args:
        `qvals` (numpy array or str)
            List of Q-values, assumed to be Phred scores.
            For how they are encoded, see `encoding`.
        `encoding` (str)
            If it is "numbers" then `qvals` should be a
            numpy array giving the Q-values. If it is
            "sanger", then `qvals` is a string, with 
            the score being the ASCI value minus 33.

    Returns:
        A number giving the average accuracy, or 
        `nan` if `qvals` is empty.

    Note that the probability :math:`p` of an error at a
    given site is related to the Q-value :math:`Q` by
    :math:`Q = -10 \log_{10} p`.

    >>> qvals = numpy.array([13, 77, 93])
    >>> round(qvalsToAccuracy(qvals), 3)
    0.983
    >>> round(qvalsToAccuracy(qvals[1 : ]), 3)
    1.0
    >>> qvalsToAccuracy(numpy.array([]))
    nan

    >>> qvals = '.n~'
    >>> round(qvalsToAccuracy(qvals, encoding='sanger'), 3)
    0.983
    """
    if len(qvals) == 0:
        return numpy.nan

    if encoding == 'numbers':
        pass
    elif encoding == 'sanger':
        qvals = numpy.array([ord(q) - 33 for q in qvals])
    else:
        raise RuntimeError("invalid `encoding`: {0}".format(encoding))

    return (1 - 10**(qvals / -10)).sum() / len(qvals)


def summarizeCCSreports(ccslist, report_type, plotfile,
                        plotminfrac=0.005):
    """Summarize and plot `CCS` reports.

    Args:
        `ccslist` (`CCS` object or list of them)
            `CCS` objects to summarize
        `report_type` (str "zmw" or "subread")
            Which type of report to summarize
        `plotfile` (str)
            Name of created bar plot
        `plotminfrac` (float)
            Only plot status categories with >=
            this fraction in at least one `CCS`

    Returns:
        Returns a pandas DataFrame aggregating the reports,
        and creates `plotfile`.
    """
    if isinstance(ccslist, CCS):
        ccslist = [ccslist]
    assert all([isinstance(ccs, CCS) for ccs in ccslist]), \
            "`ccslist` not a list of `CCS` objects"

    assert report_type in ['zmw', 'subread']
    report = report_type + '_report'

    df = (pandas.concat([getattr(ccs, report).assign(sample=ccs.samplename)
                for ccs in ccslist])
          .sort_values(['sample', 'number'], ascending=False)
          [['sample', 'status', 'number', 'fraction']]
          )

    # version of df that only has categories with `plotminfrac`
    plot_df = (df.assign(maxfrac=lambda x: x.groupby('status')
                         .fraction.transform('max'))
                 .query('maxfrac >= @plotminfrac')
                 )
    nstatus = len(plot_df.status.unique())

    p = (ggplot(plot_df) +
            geom_col(aes(x='sample', y='number', fill='status'),
                     position='stack') +
            theme(axis_text_x=element_text(angle=90, vjust=1,
                  hjust=0.5)) +
            ylab({'zmw':'ZMWs', 'subread':'subreads'}[report_type])
            )
    
    if nstatus <= len(COLOR_BLIND_PALETTE):
        p = p + scale_fill_manual(list(reversed(
                COLOR_BLIND_PALETTE[ : nstatus])))
    p.save(plotfile, 
           height=3,
           width=(2 + 0.3 * len(ccslist)),
           verbose=False)
    plt.close()

    return df

def re_expandIUPAC(re_str):
    """Expand IUPAC ambiguous nucleotide codes in `re` search string.

    Simplifies writing `re` search strings that include ambiguous
    nucleotide codes.

    Args:
        `re_str` (str)
            String appropriate to be passed to `re.compile`.

    Returns:
        A version of `re_str` where any characters not in the group
        names that correspond to upper-case ambiguous nucleotide codes
        are expanded according to their definitions in the
        `IUPAC code <https://en.wikipedia.org/wiki/Nucleic_acid_notation>`_.

    >>> re_str = '^(?P<termini5>ATG)(?P<cDNA>N+)A+(?P<barcode>N{4})$'
    >>> re_expandIUPAC(re_str)
    '^(?P<termini5>ATG)(?P<cDNA>[ACGT]+)A+(?P<barcode>[ACGT]{4})$'
    """
    # We simply do a simple replacement on all characters not in group
    # names. So first we must find group names:
    groupname_indices = set([])
    groupname_matcher = re.compile('\(\?P<[^>]*>')
    for m in groupname_matcher.finditer(re_str):
        for i in range(m.start(), m.end()):
            groupname_indices.add(i)
    
    # now replace ambiguous characters
    new_re_str = []
    for i, c in enumerate(re_str):
        if (i not in groupname_indices) and c in dms_tools2.NT_TO_REGEXP:
            new_re_str.append(dms_tools2.NT_TO_REGEXP[c])
        else:
            new_re_str.append(c)

    return ''.join(new_re_str)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
