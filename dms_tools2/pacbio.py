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

import numpy
import pandas
import pysam

from dms_tools2 import NT_TO_REGEXP

# import dms_tools2.plot to set plotting contexts / themes
import dms_tools2.plot
from dms_tools2.plot import COLOR_BLIND_PALETTE
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
        `name` (str)
            Sample or sequencing run name
        `bamfile` (str)
            BAM file created by ``ccs``
        `reportfile` (str or `None`)
            Report file created by ``ccs``, or
            `None` if you have no reports.

    Attributes:
        `name` (str)
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

    def __init__(self, name, bamfile, reportfile):
        """See main class doc string."""
        self.name = name

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


    def plotResults(self, plotfile,
            cols=['passes', 'CCS_accuracy', 'CCS_length'],
            title=True):
        """Plots the CCS results in `df`.

        The plot shows the distribution of each variable
        as well as all pairwise correlations.

        Args:
            `plotfile` (str)
                Name of created plot.
            `cols` (list)
                List of variables to plot. There must be a
                column for each in `df`.
            `title` (bool or str)
                If `False`, no title. If `True`, make
                `name` the title. If a string, make
                that the title.
        """
        assert set(cols) <= set(self.df.columns)

        def hist1d(x, color, **kwargs):
            """1D histogram for diagonal elements."""
            plt.hist(x, color=color,
                    bins=dms_tools2.plot.hist_bins_intsafe(x),
                    **kwargs)

        def hist2d(x, y, color, **kwargs):
            """2D histogram for off-diagonal elements."""
            plt.hist2d(x, y, 
                    bins=[dms_tools2.plot.hist_bins_intsafe(a)
                          for a in [x, y]],
                    cmap=dms_tools2.plot.from_white_cmap(color),
                    **kwargs)

        color = COLOR_BLIND_PALETTE[0]
        g = (dms_tools2.plot.AugmentedPairGrid(self.df, vars=cols,
                diag_sharey=False, size=3)
             .map_diag(hist1d, color=color)
             .map_upper(hist2d, color=color)
             .map_lower(hist2d, color=color)
             .ax_lims_clip_outliers()
             )

        if title:
            if not isinstance(title, str):
                title = self.name
            g.fig.suptitle(title, va='bottom')

        g.savefig(plotfile)
        plt.close()


    def filterSeqs(self, match_str, filter_colname='pass_filter',
                   expandIUPAC=True):
        """Identify CCSs that match a specific pattern.

        This filtering is useful if the CCS sequences in `df`
        should have specific subsequences.

        Args:
            `match_str` (str)
                A string that can be passed to `re.compile` that
                gives the pattern that we are looking for, with
                target subsequences as named groups. See also
                the `expandIUPAC` parameter.
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

        After calling this function, the following columns are added
        to `df`:
            Column with name given by `filter_colname`:
                Entries `True` when `match_str` matches for that row.
            Column with name given by `filter_colname` suffixed by "_polarity":
                This entry is "plus" or "minus" depending on whether we 
                directly match CCS (plus) or have to reverse-complement
                first (minus).
            Columns with name of each group in `match_str`:
                The string that matches that group in CCS, or an 
                empty string if no match. If the match is in the minus
                polarity, the string placed in this column is
                reverse-complemented to be in the plus orientation.
            Columns with group names suffixed by "_qvals":
                The Q-values for that group, or an empty numpy
                array if no match. If the group is reverse-complement
                to plus orientation, the Q-values are also reversed
                to keep them in the same orientation as the group sequence.
            Columns with group names suffixed by "_accuracy":
                The accuracy for that group, or 0 if no match.
        """
        assert filter_colname not in self.df.columns,\
                "`df` already has column {0}".format(filter_colname)
        polarity_colname = filter_colname + "_polarity"
        assert polarity_colname + "_polarity" not in self.df.columns,\
                "`df` already has column {0}".format(polarity_colname)

        if expandIUPAC:
            all_nts = ''.join(NT_TO_REGEXP.keys())
            assert not re.search('<[^>]*[ACGT]+[^>]*>', match_str),\
                    "`match_str` group name has nucleotide code"
            match_str = match_str.translate(str.maketrans(NT_TO_REGEXP))

        matcher = re.compile(match_str)
        groupnames = set(matcher.groupindex.keys())
        groupqvals = {g + '_qvals' for g in groupnames}
        groupaccuracies = {g + '_accuracy' for g in groupnames}

        # make sure created columns don't already exist
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
        assert set(match_d.keys()).isdisjoint(self.df.columns)
        self.df = pandas.concat(
                [self.df,
                 pandas.DataFrame(match_d).set_index(indexname),
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

    df = (pandas.concat([getattr(ccs, report).assign(sample=ccs.name)
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
