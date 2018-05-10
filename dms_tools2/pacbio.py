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
        `reportfile` (str)
            Report file created by ``ccs``

    Attributes:
        `name` (str)
            Name set at initialization
        `bamfile` (str)
            ``ccs`` BAM file set at initialization
        `reportfile` (str)
            ``ccs`` report file set at initialization
        `zmw_report` (pandas.DataFrame): 
            ZMW stats in `reportfile`.
            Columns are *status*, *number*, *percent*, and *fraction*.
        `subread_report` (pandas.DataFrame)
            Like `zmw_report` but for subreads.
        `df` (pandas.DataFrame)
            The CCSs in `bamfile`. Each row is a different CCS
            On creation, there will be the following columns (you
            can modify to add more): *CCS*, *CCS_qvals*, *name*,
            *passes*, *CCS_accuracy*, *CCS_length*.
    """

    def __init__(self, name, bamfile, reportfile):
        """See main class doc string."""
        self.name = name

        assert os.path.isfile(bamfile), "can't find {0}".format(bamfile)
        self.bamfile = bamfile

        assert os.path.isfile(reportfile), "can't find {0}".format(reportfile)
        self.reportfile = reportfile

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

        After calling this function, `df` has been updated by
        adding the column specified by `colname`. In addition,
        columns are added with the name of each group in 
        `match_str`. The column is an empty string if there
        is no match for the CCS in that row; otherwise it
        is the string that matches that group in CCS.
        For each group, we also add a column suffixed with 
        "_qvals" that has the Q-values, and a column suffixed
        with "_accuracy" that gives the accuracy as calculated
        from the Q-values. The "_qvals" column is an empty
        numpy array if no match, and the "_accuracy" is 0
        if no match.
        """
        assert filter_colname not in self.df.columns,\
                "`df` already has column {0}".format(filter_colname)

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
                groupaccuracies])}
        for tup in self.df.itertuples():
            s = getattr(tup, 'CCS')
            qs = getattr(tup, 'CCS_qvals')
            m = matcher.search(s)
            if m:
                match_d[filter_colname] = True
                for g in groupnames:
                    match_d[g] = m.group(g)
                    match_d[g + '_qvals'] = qs[m.start(g) : m.end(g)]
                    match_d[g + '_accuracy'] = qvalsToAccuracy(
                            match_d[g + '_qvals'])
            else:
                match_d[filter_colname] = False
                for c in groupnames:
                    match_d[c] = ''
                for c in groupqvals:
                    match_d[c] = numpy.array([], dtype='int')
                for c in groupaccuracies:
                    match_d[c] = 0

        # set index to make sure matches `df`
        indexname = df.index.name
        assert indexname not in match_d
        match_d[indexname] = df.index.tolist()
        assert set(match_d.keys()).isdisjoint(df.columns)
        self.df = pandas.concat([
                [self.df,
                 pandas.DataFrame(match_d).set_index(indexname)],
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


def qvalsToAccuracy(qvals):
    """Converts set of quality scores into average accuracy.

    Args:
        `qvals` (numpy array)
            List of Q-values, assumed Sanger encoding.

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
    """
    if len(qvals) == 0:
        return numpy.nan
    else:
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
