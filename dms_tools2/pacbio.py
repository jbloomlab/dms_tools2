"""
===================
pacbio
===================

Tools for processing PacBio sequencing data.

Some functionality of this module requires you to
have `samtools <http://www.htslib.org/>`_ installed.
"""


import os
import re
import io
import subprocess

import pandas

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

    Args:
        `name` (str): sample or sequencing run
        `bamfile` (str): ``*_ccs.bam`` file created by ``ccs``
        `reportfile` (str): ``*_report.txt`` created by ``ccs``

    Attributes:
        `name` (str)
        `bamfile` (str)
        `reportfile` (str)
        `zmw_report` (pandas.DataFrame): ZMW stats in `reportfile`.
            Columns are *status*, *number*, *percent*, and *fraction*.
        `subread_report`: Like `zmw_report` but for subreads.
    """

    def __init__(self, name, bamfile, reportfile):
        """See main class doc string."""
        self.name = name

        assert os.path.isfile(bamfile), "can't find {0}".format(bamfile)
        self.bamfile = bamfile

        assert os.path.isfile(reportfile), "can't find {0}".format(reportfile)
        self.reportfile = reportfile

        self._samfile = None

        # set `zmw_report` and `subread_report`
        self._parse_report()


    @property
    def samfile(self):
        """Create and return SAM version of `bamfile`."""
        if (self._samfile is None) or (not os.path.isfile(self._samfile)):
            try:
                subprocess.check_output(['samtools', '--version'])
            except:
                raise RuntimeError("`samtools` must be installed")
            self._samfile = os.path.splitext(self.bamfile)[0] + '.sam'
            try:
                subprocess.check_call(['samtools', 'view', 
                        '-o', self._samfile, self.bamfile])
            except:
                if os.path.isfile(samfile):
                    os.remove(samfile) # avoid incomplete files
                raise
            assert os.path.isfile(self._samfile), \
                    "Failed to create {0}".format(self._samfile)

        return self._samfile


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



def summarizeCCSreports(ccslist, report_type, plotfile,
                        plotminfrac=0.005):
    """Summarize and plot `CCS` reports.

    Args:
        `ccslist`: a `CCS` object or a list of them
        `report_type` (str): "zmw" or "subread" indicating
            which type of report to summarize
        `plotfile` (str): name of created bar plot
        `plotminfrac` (float): only plot status categories
            wth >= this fraction in at least one `CCS`

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
           width=(2.5 + 0.3 * len(ccslist)),
           verbose=False)

    return df


if __name__ == '__main__':
    import doctest
    doctest.testmod()
