"""
===================
plot
===================

Plotting functions for ``dms_tools2``.
"""


import os
import pandas
from plotnine import *
theme_set(theme_classic()) # classic ggplot theme


def plotReadStats(names, readstatfiles, plotfile):
    """Plots `dms2_bcsubamp` read statistics for a set of samples.
    
    Args:
        `names` (list or series)
            Names of the samples for which we are plotting statistics.
        `readstatfiles` (list or series)
            Names of ``*_readstats.csv`` files created by ``dms2_bcsubamp``
            containing the data to plot.
        `plotfile` (str)
            Name of PDF plot file to create.
    """
    assert len(names) == len(readstatfiles)
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'
    readstats = pandas.concat([pandas.read_csv(f).assign(name=name) for
                (name, f) in zip(names, readstatfiles)], ignore_index=True)
    readstats['retained'] = (readstats['total'] - readstats['fail filter']
            - readstats['low Q barcode'])
    readstats_melt = readstats.melt(id_vars='name', 
            value_vars=['retained', 'fail filter', 'low Q barcode'],
            value_name='number of reads', var_name='read fate')
    p = (ggplot(readstats_melt)
            + geom_col(aes(x='name', y='number of reads', fill='read fate'),
                position='stack')
            + theme(axis_text_x=element_text(angle=90, vjust=1, hjust=0.5),
                    axis_title_x=element_blank()) 
            )
    p.save(plotfile, height=2.7, width=(1.2 + 0.4 * len(names)))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
