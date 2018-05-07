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
            can modify to add more): *CCS*, *qvals*, *name*,
            *passes*, *accuracy*, *length*.
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
            cols=['passes', 'accuracy', 'length'],
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
            d['qvals'].append(s.query_qualities)
            d['name'].append(s.query_name)
            d['passes'].append(s.get_tag('np'))
            d['accuracy'].append(s.get_tag('rq'))
            d['length'].append(s.query_length)

        # create data frame
        self.df = pandas.DataFrame(d)

        # some checks on `df`
        assert self.df.name.size == self.df.name.unique().size,\
                "non-unique names for {0}".format(self.name)
        assert (self.df.length == self.df.CCS.apply(len)).all(),\
                "CCS not correct length"
        assert (self.df.length == self.df.qvals.apply(len)).all(),\
                "qvals not correct length"



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
