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

from dms_tools2 import NT_TO_REGEXP

# import dms_tools2.plot to set plotting contexts / themes
import dms_tools2.plot
from dms_tools2.plot import COLOR_BLIND_PALETTE
from dms_tools2.plot import COLOR_BLIND_PALETTE_GRAY
import matplotlib.pyplot as plt
import seaborn
from plotnine import *


class CCS:
    """Class to handle results of ``ccs``.

    Reads and manipulate results of PacBio ``ccs``.
    Has been tested on output of ``ccs`` version 3.0.0.

    This class reads all data into memory, and so you
    may need a lot of RAM if `bamfile` is large.

    Args:
        `sample` (str)
            Sample or sequencing run
        `bamfile` (str)
            BAM file created by ``ccs``
        `reportfile` (str or `None`)
            Report file created by ``ccs``, or
            `None` if you have no reports.

    Attributes:
        `sample` (str)
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
            can modify to add more): *CCS*, *CCS_qvals*, *name*,
            *passes*, *CCS_accuracy*, *CCS_length*.

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

    Check `ccs.df` has correct names, CCS sequences,
    and columns:

    >>> set(ccs.df.name) == {s['name'] for s in ccs_seqs}
    True
    >>> set(ccs.df.CCS) == {s['seq'] for s in ccs_seqs}
    True
    >>> set(ccs.df.columns) == {'CCS', 'CCS_qvals', 'name',
    ...         'passes', 'CCS_accuracy', 'CCS_length'}
    True

    Apply filter for sequences that have expected termini
    and define barcode and read in these:

    >>> match_str = (termini5 + '(?P<barcode>N{3})' +
    ...         '(?P<read>N+)' + termini3)
    >>> ccs.filterSeqs(match_str, filter_colname='barcoded')

    This filtering add new columns to `ccs.df`:

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
    >>> numpy.allclose(ccs.df.query('not barcoded').barcode_accuracy,
    ...     0, atol=1e-4)
    True
    >>> barcoded_polarity = [x.split('_')[3] for x in bc_names]
    >>> (ccs.df.query('barcoded').barcoded_polarity == barcoded_polarity).all()
    True

    """

    def __init__(self, sample, bamfile, reportfile):
        """See main class doc string."""
        self.sample = sample

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


    def plotResults(self, plotfile, lower_filter=None,
            cols=['passes', 'CCS_accuracy', 'CCS_length'],
            title=True):
        """Plots the CCS results in `df`.

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
                `sample` the title. If a string, make
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
                title = self.sample
            g.fig.suptitle(title, va='bottom')

        g.savefig(plotfile)
        plt.close()


    def filterSeqs(self, match_str, filter_colname='pass_filter',
                   expandIUPAC=True, overwrite=True):
        """Identify CCSs that match a specific pattern.

        This filtering is useful if the CCS sequences in `df`
        should have specific subsequences.

        Args:
            `match_str` (str)
                A string that can be passed to `re.compile` that
                gives the pattern that we are looking for, with
                target subsequences as named groups. See also
                the `expandIUPAC` parameter. Note that if 
                `expandIUPAC` is true, the group names cannot
                include IUPAC codes (a safe thing to do is just
                make the group names all lower case).
            `filter_colname` (str)
                Name of a new column added to `df`. Every row
                in `df` that has a CCS column that matches
                `match_str` gets this column set to `True`,
                and all other rows get it set to `False`.
            `expandIUPAC` (bool)
                Use `IUPAC code <https://en.wikipedia.org/wiki/Nucleic_acid_notation>`_
                to expand ambiguous nucleotides (e.g., "N") in
                `match_str`. This can simplify the writing of
                `match_str`. If you use this option, ensure that
                none of the named groups in `match_str` have upper-
                case letters that are nucleotide codes.
            `overwrite` (bool)
                If `True`, we overwrite any existing columns to
                be created that already exist. If `False`, raise
                an error if any of the columns already exist.

        The following columns are added to `df`:

            - Name given by `filter_colname`: `True` if `match_str` matches.

            - Columns with name of `filter_colname` and these suffixes:
                - "_polarity": "plus" or "minus" depending on whether 
                  matches CCS directly or matches reverse-complement.
                - Each group in `match_str`: the string that matches
                  that group in CCS, or any empty string if no match.
                  If the match is in the minus polarity, the string is
                  reverse-complemented to be in the plus orientation.
                - Each group names suffixed by "_qvals": The Q-values
                  for that group, or an empty numpy if no match. If
                  the group is reverse-complemented to plus orientation,
                  the Q-values are also reversed to keep them in the
                  same orientation as the group sequence.
                - Group names suffixed by "_accuracy": accuracy for
                  that group, or 0 if no match.
        """
        polarity_colname = filter_colname + "_polarity"

        if expandIUPAC:
            expand_nts = ''.join([key for key, value in 
                    NT_TO_REGEXP.items() if len(value) > 1])
            assert not re.search('<[^>]*[{0}]+[^>]*>'.format(expand_nts),
                    match_str), "`match_str` group name has IUPAC code"
            match_str = match_str.translate(str.maketrans(NT_TO_REGEXP))

        matcher = re.compile(match_str)
        groupnames = set(matcher.groupindex.keys())
        groupqvals = {g + '_qvals' for g in groupnames}
        groupaccuracies = {g + '_accuracy' for g in groupnames}

        # make sure created columns don't already exist
        if not overwrite:
            assert filter_colname not in self.df.columns,\
                    "`df` already has column {0}".format(filter_colname)
            assert polarity_colname + "_polarity" not in self.df.columns,\
                    "`df` already has column {0}".format(polarity_colname)
            assert groupnames.isdisjoint(self.df.columns), \
                    "`df` has columns with `match_str` group names"
            assert groupqvals.isdisjoint(self.df.columns), \
                    "`df` has columns with `match_str` group qvals"
            assert groupaccuracies.isdisjoint(self.df.columns), \
                    "`df` has columns with `match_str` group accuracies"

        # look for matches for each row
        match_d = {c:[] for c in set.union(*[groupnames, groupqvals,
                groupaccuracies, {filter_colname, polarity_colname}])}
        for tup in self.df.itertuples():
            s = getattr(tup, 'CCS')
            qs = getattr(tup, 'CCS_qvals')
            m = matcher.search(s)
            if m:
                polarity = "plus"
            else:
                m = matcher.search(dms_tools2.utils.reverseComplement(s))
                qs = numpy.flip(qs, axis=0)
                polarity = "minus"
            if m:
                match_d[filter_colname].append(True)
                match_d[polarity_colname].append(polarity)
                for g in groupnames:
                    if polarity == "plus":
                        match_d[g].append(m.group(g))
                    else:
                        assert polarity == "minus"
                        match_d[g].append(m.group(g))
                    g_qvals = qs[m.start(g) : m.end(g)]
                    match_d[g + '_qvals'].append(g_qvals)
                    match_d[g + '_accuracy'].append(
                            qvalsToAccuracy(g_qvals))
            else:
                match_d[filter_colname].append(False)
                match_d[polarity_colname].append('')
                for c in groupnames:
                    match_d[c].append('')
                for c in groupqvals:
                    match_d[c].append(numpy.array([], dtype='int'))
                for c in groupaccuracies:
                    match_d[c].append(0)

        # set index to make sure matches `df`
        indexname = self.df.index.name
        assert indexname not in match_d
        match_d[indexname] = self.df.index.tolist()
        dup_cols = set(match_d.keys()).intersection(set(self.df.columns))
        if (not overwrite) and dup_cols:
            raise ValueError("overwriting columns")
        self.df = pandas.concat(
                [self.df.drop(dup_cols, axis=1),
                 pandas.DataFrame(match_d).set_index(indexname),
                ],
                axis=1)

    def align(self, mapper, query_col, alignment_col='aligned',
              overwrite=True, paf_file=None):
        """Align sequences to target sequence(s).

        Arguments:
            `mapper` (:py:mod:`dms_tools2.minimap2.Mapper`)
                Align sequences using `mapper.map`. Target
                sequences are specified when initializing `mapper`.
            `query_col` (str)
                Column in `df` with query sequences to align.
            `alignment_col` (str)
				Specify names of new columns added to `df` (see
				below).
            `overwrite` (bool)
                If `True`, we overwrite any existing columns to
                be created that already exist. If `False`, raise
                an error if any of the columns already exist.
            `paf_file` (`None` or str)
                If a str, is the name of the PAF file created
                by `mapper` (see `outfile` argument to
                :py:mod:`dms_tools2.minimap2.Mapper.map`).

        Calling this function adds the following columns to `df`:

            - Name given by `alignment_col`: `True` if alignment.

            - Columns with names of `alignment_col` and these suffixes:
                - "_target": name of target to which query aligned,
                  or empty string if no alignment.
                - "_clip_start": number of nucleotides clipped from
                  start of query, or -1 if no alignment.
                - "_clip_end": number of nucleotides clipped from
                  end of query, or -1 if no alignment.
                - "_start": position in target where alignment starts in
                  0-based indexing, or -1 if no alignment.
                - "_cigar" Long format cigar string 
                  (`see here <https://github.com/lh3/minimap2>`_), 
                  or empty string if no alignment
        """
        assert query_col in self.df.columns, \
                "no `query_col` {0}".format(query_col)

        newcols = {'aligned':alignment_col,
                   'target':alignment_col + '_target',
                   'clip_start':alignment_col + '_clip_start',
                   'clip_end':alignment_col + '_clip_end',
                   'start':alignment_col + '_start',
                   'cigar':alignment_col + '_cigar'
                  }
        dup_cols = set(newcols.keys()).intersection(set(self.df.columns))
        if (not overwrite) and dup_cols:
            raise ValueError("would duplicate existing columns:\n{0}"
                             .format(dup_cols))

        assert len(self.df.name) == len(self.df.name.unique()), \
                "`name` in `df` not unique"
        with tempfile.NamedTemporaryFile(mode='w') as queryfile:
            queryfile.write('\n'.join([
                            '>{0}\n{1}'.format(*tup) for tup in
                            self.df.query('{0} != ""'.format(query_col))
                                [['name', query_col]]
                                .itertuples(index=False, name=False)
                            ]))
            map_dict = mapper.map(queryfile.name, outfile=paf_file)

        align_d = {c:[] for c in newcols.values()}
        for name in self.df.name:
            if name in map_dict:
                a = map_dict[name]
                assert a.strand == 1, "method does not handle - polarity"
                align_d[newcols['aligned']].append(True)
                align_d[newcols['target']].append(a.target)
                align_d[newcols['cigar']].append(a.cigar_str)
                align_d[newcols['start']].append(a.r_st)
                align_d[newcols['clip_start']].append(a.q_st)
                align_d[newcols['clip_end']].append(a.q_len - a.q_en)

            else:
                align_d[newcols['aligned']].append(False)
                for col in ['target', 'cigar']:
                    align_d[newcols[col]].append('')
                for col in ['clip_start', 'clip_end', 'start']:
                    align_d[newcols[col]].append('')

        # add to df, making sure index is correct
        index_name = self.df.index.name
        assert index_name not in align_d
        align_d[index_name] = self.df.index.tolist()
        if (not overwrite) and dup_cols:
            raise ValueError("overwriting columns")
        self.df = pandas.concat(
                [self.df.drop(dup_cols, axis=1),
                 pandas.DataFrame(align_d).set_index(index_name),
                ],
                axis=1)


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

        # create data frame
        self.df = pandas.DataFrame(d)

        # some checks on `df`
        assert self.df.name.size == self.df.name.unique().size,\
                "non-unique names for {0}".format(self.name)
        assert (self.df.CCS_length == self.df.CCS.apply(len)).all(),\
                "CCS not correct length"
        assert (self.df.CCS_length == self.df.CCS_qvals.apply(len)).all(),\
                "qvals not correct length"


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

    df = (pandas.concat([getattr(ccs, report).assign(sample=ccs.sample)
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


if __name__ == '__main__':
    import doctest
    doctest.testmod()
