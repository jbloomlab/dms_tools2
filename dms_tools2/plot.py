"""
===================
plot
===================

Plotting functions for ``dms_tools2``.

The plotting is done with ``plotnine``.

The default ``plotnine`` color palette is **not** color-blind safe.
Therefore, this module defines the constant `COLOR_BLIND_PALETTE`
which gives the color-blind safe scale described here:
http://bconnelly.net/2013/10/creating-colorblind-friendly-figures/

You can use this scale by adding the following to your plots:
`scale_fill_manual(COLOR_BLIND_PALETTE)` or
`scale_color_manual(COLOR_BLIND_PALETTE)`.
"""


import os
import math
import pandas
from plotnine import *
theme_set(theme_classic()) # use classic ggplot theme
import dms_tools2.utils

COLOR_BLIND_PALETTE = ["#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]


def latexSciNot(xlist, exp_cutoff=3, ddigits=1):
    """Converts list of numbers to LaTex scientific notation.

    Useful for nice axis-tick formatting.

    Args:
        `xlist` (list)
            Numbers to format.
        `exp_cutoff` (int)
            Convert to scientific notation if `abs(math.log10(x))` >= this.
        `ddigits` (int)
            Show at most this many digits after the decimal place, shows
            less if not needed to precisely express all numbers.

    Returns:
        List of latex scientific notation formatted strings.

    >>> latexSciNot([0, 3, 3120, 0.07, 0.000927])
    ['$0$', '$3.0$', '$3.1 \\\\times 10^{3}$', '$0.1$', '$9.3 \\\\times 10^{-4}$']

    >>> latexSciNot([0.001, 1, 1000, 1e6])
    ['$10^{-3}$', '$1$', '$10^{3}$', '$10^{6}$']
    """
    # can all numbers be expressed as 10**integer?
    if all([10**int(math.log10(x)) == x for x in xlist if x != 0]):
        all_exp10 = True
    else:
        all_exp10 = False

    # can all numbers be expressed as integer * 10**integer?
    if all([x == (10**int(math.log10(x))) * int(x / 10**(int(math.log10(x))))
            for x in xlist if x != 0]):
        ddigits = 0

    # make formatted numbers
    formatlist = []
    for x in xlist:
        if x < 0:
            raise ValueError("only handles numbers >= 0")
        elif x == 0:
            formatlist.append('$0$')
            continue
        exponent = int(math.log10(x))
        if math.log10(x) < exponent and x < 1:
            exponent -= 1
        if all_exp10:
            if abs(exponent) >= exp_cutoff:
                xformat = '10^{{{0}}}'.format(exponent)
            else:
                xformat = str(int(x))
        elif abs(exponent) >= exp_cutoff:
            formatstr = '{0:.' + str(ddigits) + 'f} \\times 10^{{{1}}}'
            xformat = formatstr.format(x / 10.**exponent, exponent)
        else:
            formatstr = '{0:.' + str(ddigits) + 'f}'
            xformat = formatstr.format(x)
        formatlist.append('${0}$'.format(xformat))
    return formatlist


def plotReadStats(names, readstatfiles, plotfile):
    """Plots `dms2_bcsubamp` read statistics for a set of samples.
    
    Args:
        `names` (list or series)
            Names of the samples for which we are plotting statistics.
        `readstatfiles` (list or series)
            Names of ``*_readstats.csv`` files created by ``dms2_bcsubamp``.
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
            + scale_y_continuous(labels=latexSciNot) 
            + scale_fill_manual(COLOR_BLIND_PALETTE)
            )
    p.save(plotfile, height=2.7, width=(1.2 + 0.3 * len(names)))


def plotBCStats(names, bcstatsfiles, plotfile):
    """Plots `dms2_bcsubamp` barcode statistics for a set of samples.

    Args:
        `names` (list or series)
            Names of the samples for which we are plotting statistics.
        `bcstatsfiles` (list or series)
            Names of ``*_bcstats.csv`` files created by ``dms2_bcsubamp``.
        `plotfile` (str)
            Name of PDF plot file to create.
    """
    assert len(names) == len(bcstatsfiles)
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'
    bcstats = pandas.concat([pandas.read_csv(f).assign(name=name) for
                (name, f) in zip(names, bcstatsfiles)], ignore_index=True)
    bcstats_melt = bcstats.melt(id_vars='name', 
            value_vars=['too few reads', 'not alignable', 'aligned'],
            value_name='number of barcodes', var_name='barcode fate')
    p = (ggplot(bcstats_melt)
            + geom_col(aes(x='name', y='number of barcodes', 
                fill='barcode fate'), position='stack')
            + theme(axis_text_x=element_text(angle=90, vjust=1, hjust=0.5),
                    axis_title_x=element_blank())
            + scale_y_continuous(labels=latexSciNot)
            + scale_fill_manual(COLOR_BLIND_PALETTE)
            )
    p.save(plotfile, height=2.7, width=(1.2 + 0.3 * len(names)))


def plotReadsPerBC(names, readsperbcfiles, plotfile, 
        maxreads=10, maxcol=6):
    """Plots `dms2_bcsubamp` reads-per-barcode stats for set of samples.

    Args:
        `names` (list or series)
            Names of samples for which we plot statistics.
        `readsperbcfiles` (list or series)
            Names of ``*_readsperbc.csv`` files created by ``dms2_bcsubamp``.
        `plotfile` (str)
            Name of PDF plot file to create.
        `maxreads` (int)
            For any barcodes with > this many reads, just make a category
            of >= this.
        `maxcol` (int)
            Number of columns in faceted plot.
    """
    assert len(names) == len(readsperbcfiles)
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'

    # read data frames, ensure 'number of reads' from 1 to >= maxreads
    dfs = []
    for (name, f) in zip(names, readsperbcfiles):
        df = pandas.read_csv(f)
        # make 'number of reads' maxreads hold number >= maxreads barcodes
        n_ge = df[df['number of reads'] >= maxreads]['number of reads'].sum()
        df = df.append(pandas.DataFrame({'number of reads':[maxreads],
                'number of barcodes':[n_ge]}))
        for nreads in range(1, maxreads):
            if nreads not in df['number of reads']:
                df.append(pandas.DataFrame({'number of reads':[nreads],
                        'number of barcodes':[0]}))
        df = df[df['number of reads'] <= maxreads]
        df = df.assign(name=name)
        dfs.append(df)
    df = pandas.concat(dfs, ignore_index=True)

    ncol = min(maxcol, len(names))
    nrow = int(math.ceil(len(names) / float(ncol)))
    p = (ggplot(df)
            + geom_col(aes(x='number of reads', y='number of barcodes'),
                position='stack')
            + scale_x_continuous(breaks=[1, maxreads // 2, maxreads],
                    labels=['$1$', '${0}$'.format(maxreads // 2), 
                    '$\geq {0}$'.format(maxreads)])
            + scale_y_continuous(labels=latexSciNot)
            + facet_wrap('~name', ncol=ncol) 
            + theme(strip_text=element_text(lineheight=1.8)))

    p.save(plotfile, 
            height=1.2 * (0.4 + nrow),
            width=(1.5 * (0.8 + ncol)))


def plotDepth(names, countsfiles, plotfile, maxcol=5):
    """Plot sequencing depth along primary sequence.

    Args:
        `names` (list or series)
            Names of samples for which we plot statistics.
        `countsfiles` (list or series)
            ``*_codoncounts.csv`` files of type created by ``dms2_bcsubamp``.
        `plotfile` (str)
            Name of created PDF plot file containing count depth.
        `maxcol` (int)
            Number of columns in faceted plot.
    """
    assert len(names) == len(countsfiles)
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'

    counts = pandas.concat([dms_tools2.utils.annotateCodonCounts(f).assign(
            name=name) for (name, f) in zip(names, countsfiles)], 
            ignore_index=True)
    
    ncol = min(maxcol, len(names))
    nrow = int(math.ceil(len(names) / float(ncol)))

    p = (ggplot(counts, aes(x='site', y='ncounts'))
            + geom_line()
            + scale_y_continuous(labels=latexSciNot, 
                    limits=(0, counts['ncounts'].max()), name='number of counts')
            + scale_x_continuous(limits=(counts['site'].min(),
                    counts['site'].max()))
            + facet_wrap('~name', ncol=ncol) 
            + theme(strip_text=element_text(lineheight=1.8)))

    p.save(plotfile, 
            height=1.2 * (0.4 + nrow),
            width=(1.8 * (0.6 + ncol)))


def plotMutFreq(names, countsfiles, plotfile, maxcol=5):
    """Plot mutation frequency along primary sequence.

    Args:
        `names` (list or series)
            Names of samples for which we plot statistics.
        `countsfiles` (list or series)
            ``*_codoncounts.csv`` files of type created by ``dms2_bcsubamp``.
        `plotfile` (str)
            Name of created PDF plot file.
        `maxcol` (int)
            Number of columns in faceted plot.
    """
    assert len(names) == len(countsfiles)
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'

    counts = pandas.concat([dms_tools2.utils.annotateCodonCounts(f).assign(
            name=name) for (name, f) in zip(names, countsfiles)], 
            ignore_index=True)
    
    ncol = min(maxcol, len(names))
    nrow = int(math.ceil(len(names) / float(ncol)))

    p = (ggplot(counts, aes(x='site', y='mutfreq'))
            + geom_line()
            + scale_y_continuous(labels=latexSciNot, 
                    limits=(0, counts['mutfreq'].max()), 
                    name='mutation frequency')
            + scale_x_continuous(limits=(counts['site'].min(),
                    counts['site'].max()))
            + facet_wrap('~name', ncol=ncol) 
            + theme(strip_text=element_text(lineheight=1.8)))

    p.save(plotfile, 
            height=1.2 * (0.4 + nrow),
            width=(1.8 * (0.6 + ncol)))


def plotCodonMutTypes(names, countsfiles, plotfile,
        classification='aachange'):
    """Plot average frequency codon mutation types.

    The averages are determined by summing counts for all sites.

    Args:
        `names` (list or series)
            Names of samples for which we plot statistics.
        `countsfiles` (list or series)
            ``*_codoncounts.csv`` files of type created by ``dms2_bcsubamp``.
        `plotfile` (str)
            Name of created PDF plot file.
        `classification` (str)
            The method used to classify the mutation types. Can be:

                `aachange` : classify as stop, synonymous, nonsynonymous

                `n_ntchanges` : classify by number of nucleotide changes
    """
    assert len(names) == len(countsfiles)
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'

    counts = pandas.concat([dms_tools2.utils.annotateCodonCounts(f).assign(
            name=name) for (name, f) in zip(names, countsfiles)], 
            ignore_index=True)

    if classification == 'aachange':
        muttypes = {'stop':'nstop', 'synonymous':'nsyn', 
                'nonsynonymous':'nnonsyn'}
    elif classification == 'n_ntchanges':
        muttypes = dict([('{0} nucleotide'.format(n), 'n{0}nt'.format(n))
                for n in [1, 2, 3]])
    else:
        raise ValueError("Invalid classification {0}".format(classification))

    df = (counts[list(muttypes.values()) + ['ncounts', 'name']]
            .groupby('name', as_index=False)
            .sum(axis=1)
            .assign(ncounts=lambda x: x['ncounts'].astype('float'))
            )
    for (newcol, n) in muttypes.items():
        df[newcol] = df[n] / df['ncounts']
    df = df.melt(id_vars='name', var_name='mutation type',
                value_vars=list(muttypes.keys()),
                value_name='per-codon frequency')

    p = (ggplot(df)
            + geom_col(aes(x='name', y='per-codon frequency', 
                fill='mutation type'), position='stack')
            + theme(axis_text_x=element_text(angle=90, vjust=1, hjust=0.5),
                    axis_title_x=element_blank())
            + scale_y_continuous(labels=latexSciNot)
            + scale_fill_manual(COLOR_BLIND_PALETTE)
            )

    p.save(plotfile, height=2.7, width=(1.2 + 0.3 * len(names)))



if __name__ == '__main__':
    import doctest
    doctest.testmod()
