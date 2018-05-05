"""
===================
pacbio
===================

Tools for processing PacBio sequencing data.
"""


import os
import re
import io

import pandas

# import dms_tools2.plot to set plotting contexts / themes
import dms_tools2.plot
from dms_tools2.plot import COLOR_BLIND_PALETTE
import matplotlib.pyplot as plt
import seaborn
from plotnine import *


def summarizeCCSreport(reportfiles, samples, zmw_plotfile,
                       subread_plotfile, plotminfrac=0.005):
    """Summarize reports from ``ccs``.

    Processes ``*_report.csv`` files. Tested with formats
    from ``ccs`` version 3.0.0.

    Args:
        `reportfiles` (list or str)
            List of report files, can be string if just one.
        `samples` (list or str or `None`)
            Sample names for each report file. If `None`
            then use base file name.
        `zmw_plotfile` (str or `None`)
            Name of generated plot for ZMWs, or `None` if
            no plot should be made.
        `subread_plotfile` (str or `None`)
            Like `zmw_plotfile` but for subreads.
        `plotminfrac` (float)
            Only plot categories that have a fraction
            >= this number.

    Returns:
        A pandas DataFrame summarizing the results.
    """
    if isinstance(reportfiles, str):
        reportfiles = [reportfiles]
    if isinstance(samples, str):
        samples = [samples]
    if samples is None:
        samples = [os.path.basename(f) for f in reportfiles]
    assert len(reportfiles) == len(samples)
    assert len(samples) == len(set(samples)), "samples not unique"

    # match reports made by ccs 3.0.0
    reportmatch = re.compile('^ZMW Yield\n(?P<ZMW>(.+\n)+)\n\n'
                             'Subread Yield\n(?P<Subread>(.+\n)+)$')

    df_list = []
    for (sample, fname) in zip(samples, reportfiles):
        with open(fname) as f:
            report = f.read()
        m = reportmatch.search(report)
        assert m, "Cannot match {0}\n\n{1}".format(fname, report)
        for read_type in ['ZMW', 'Subread']:
            df_list.append(pandas.read_csv(
                                io.StringIO(m.group(read_type)),
                                names=['status', 'number', 'percent']
                                )
                           .assign(sample=sample,
                                   read_type=read_type)
                           )
    df = (pandas.concat(df_list)
          [['read_type', 'sample', 'status', 'number', 'percent']]
          .sort_values('read_type', ascending=False)
          .reset_index(drop=True)
          .assign(fraction=lambda x: x.percent
                                      .str.slice(None, -1)
                                      .astype('float')
                                      / 100)
          )

    for (plotfile, read_type) in [(zmw_plotfile, 'ZMW'), 
            (subread_plotfile, 'Subread')]:
        if plotfile is not None:
            rt_df = (df.query('read_type == @read_type')
                       .assign(maxfrac=lambda x: 
                                x.groupby('status')
                                 .fraction
                                 .transform('max'))
                       .query('maxfrac >= @plotminfrac')
                       )
            nstatus = len(rt_df.status.unique())
            p = (ggplot(rt_df) +
                    geom_col(aes(x='sample', y='number',
                        fill='status'), position='stack') +
                    theme(axis_text_x=element_text(angle=90,
                        vjust=1, hjust=0.5)) +
                    ylab(read_type + 's') 
                    )
            if nstatus <= len(COLOR_BLIND_PALETTE):
                p = p + scale_fill_manual(list(reversed(
                        COLOR_BLIND_PALETTE[ : nstatus])))
            p.save(plotfile, 
                    height=3,
                    width=(2.5 + 0.3 * len(samples)),
                    verbose=False)

    return df



if __name__ == '__main__':
    import doctest
    doctest.testmod()
