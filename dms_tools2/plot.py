"""
===================
plot
===================

Plotting related functions for ``dms_tools2``.

Uses `plotnine <https://plotnine.readthedocs.io/en/stable>`_
and `seaborn <https://seaborn.pydata.org/index.html>`_.
"""


import re
import os
import math
import numbers
import random
import collections

import natsort
import pandas
import numpy
import scipy.stats
import scipy.optimize
from statsmodels.sandbox.stats.multicomp import multipletests

# complicated backend setting: we typically want PDF,
# but this causes problem when loading from iPython / Jupyter
import matplotlib
backend = matplotlib.get_backend()
try:
    matplotlib.use('pdf', warn=False)
    import matplotlib.pyplot as plt
except:
    matplotlib.use(backend, warn=False, force=True)
    import matplotlib.pyplot as plt

from plotnine import *
# set ggplot theme
theme_set(theme_bw(base_size=12)) 

import seaborn
seaborn.set(context='talk',
            style='white',
            rc={
                'xtick.labelsize':15,
                'ytick.labelsize':15,
                'axes.labelsize':19,
                'legend.fontsize':17,
                'font.family':'sans-serif',
                'font.sans-serif':['DejaVu Sans'],
                }
           )

from dms_tools2 import CODONS, AAS, AAS_WITHSTOP, NTS
import dms_tools2.utils

#: `color-blind safe palette <http://bconnelly.net/2013/10/creating-colorblind-friendly-figures/>`_
#: use by adding to your plots the following
#: `scale_fill_manual(COLOR_BLIND_PALETTE)` or
#: `scale_color_manual(COLOR_BLIND_PALETTE)`.
COLOR_BLIND_PALETTE = ["#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

#: `color-blind safe palette <http://bconnelly.net/2013/10/creating-colorblind-friendly-figures/>`_
#: that differs from `COLOR_BLIND_PALETTE` in that first 
#: color is gray rather than black.
COLOR_BLIND_PALETTE_GRAY = ["#999999", "#E69F00", "#56B4E9", "#009E73",
                            "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]


def breaksAndLabels(xi, x, n):
    """Get breaks and labels for an axis.

    Useful when you would like to re-label a numeric x-axis
    with string labels.

    Uses `matplotlib.ticker.MaxNLocator` to choose pretty breaks.

    Args:
        `xi` (list or array)
            Integer values actually assigned to axis points.
        `x` (list)
            Strings corresponding to each numeric value in `xi`.
        `n` (int)
            Approximate number of ticks to use.

    Returns:
        The tuple `(breaks, labels)` where `breaks` gives the
        locations of breaks taken from `xi`, and `labels` is
        the label for each break.

    >>> xi = list(range(213))
    >>> x = [str(i + 1) for i in xi]
    >>> (breaks, labels) = breaksAndLabels(xi, x, 5)
    >>> breaks
    [0, 50, 100, 150, 200]
    >>> labels
    ['1', '51', '101', '151', '201']
    """
    assert len(xi) == len(x)
    assert all([isinstance(i, (int, numpy.integer)) for i in xi]), \
            "xi not integer values:\n{0}".format(xi)
    xi = list(xi)
    assert sorted(set(xi)) == xi, "xi not unique and ordered"
    breaks = matplotlib.ticker.MaxNLocator(n).tick_values(xi[0], xi[-1])
    breaks = [int(i) for i in breaks if xi[0] <= i <= xi[-1]]
    labels = [x[xi.index(i)] for i in breaks]
    return (breaks, labels)


def latexSciNot(xlist):
    """Converts list of numbers to LaTex scientific notation.

    Useful for nice axis-tick formatting.

    Args:
        `xlist` (list or number)
            Numbers to format.

    Returns:
        List of latex scientific notation formatted strings,
        or single string if `xlist` is a number.

    >>> latexSciNot([0, 3, 3120, -0.0000927])
    ['$0$', '$3$', '$3.1 \\\\times 10^{3}$', '$-9.3 \\\\times 10^{-5}$']

    >>> latexSciNot([0.001, 1, 1000, 1e6])
    ['$0.001$', '$1$', '$10^{3}$', '$10^{6}$']

    >>> latexSciNot([-0.002, 0.003, 0.000011])
    ['$-0.002$', '$0.003$', '$1.1 \\\\times 10^{-5}$']

    >>> latexSciNot([-0.1, 0.0, 0.1, 0.2])
    ['$-0.1$', '$0$', '$0.1$', '$0.2$']

    >>> latexSciNot([0, 1, 2])
    ['$0$', '$1$', '$2$']
    """
    if isinstance(xlist, numbers.Number):
        isnum = True
        xlist = [xlist]
    else:
        isnum = False
    formatlist = []
    for x in xlist:
        xf = "{0:.2g}".format(x)
        if xf[ : 2] == '1e':
            xf = "$10^{{{0}}}$".format(int(xf[2 : ]))
        elif xf[ : 3] == '-1e':
            xf = "$-10^{{{0}}}$".format(int(xf[3 : ]))
        elif 'e' in xf:
            (d, exp) = xf.split('e')
            xf = '${0} \\times 10^{{{1}}}$'.format(d, int(exp))
        else:
            xf = '${0}$'.format(xf)
        formatlist.append(xf)
    if isnum:
        assert len(formatlist) == 1
        formatlist = formatlist[0]
    return formatlist


def plotReadStats(names, readstatfiles, plotfile):
    """Plots ``dms2_bcsubamp`` read statistics for a set of samples.
    
    Args:
        `names` (list or series)
            Names of the samples for which we are plotting statistics.
        `readstatfiles` (list or series)
            Names of ``*_readstats.csv`` files created by ``dms2_bcsubamp``.
        `plotfile` (str)
            Name of PDF plot file to create.
    """
    assert len(names) == len(readstatfiles) == len(set(names))
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
    p.save(plotfile, height=2.7, width=(1.2 + 0.25 * len(names)),
            verbose=False)
    plt.close()


def plotBCStats(names, bcstatsfiles, plotfile):
    """Plots ``dms2_bcsubamp`` barcode statistics for set of samples.

    Args:
        `names` (list or series)
            Names of the samples for which we are plotting statistics.
        `bcstatsfiles` (list or series)
            Names of ``*_bcstats.csv`` files created by ``dms2_bcsubamp``.
        `plotfile` (str)
            Name of PDF plot file to create.
    """
    assert len(names) == len(bcstatsfiles) == len(set(names))
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'
    bcstats = pandas.concat([pandas.read_csv(f).assign(name=name) for
                (name, f) in zip(names, bcstatsfiles)], ignore_index=True)
    bcstats_melt = bcstats.melt(id_vars='name', 
            value_vars=['too few reads', 'not alignable', 'aligned'],
            value_name='number of barcodes', var_name='barcode fate')
    p = (ggplot(bcstats_melt)
            + geom_col(aes(x='name', y='number of barcodes', 
                fill='barcode fate'), position=position_stack(reverse=True))
            + theme(axis_text_x=element_text(angle=90, vjust=1, hjust=0.5),
                    axis_title_x=element_blank())
            + scale_y_continuous(labels=latexSciNot)
            + scale_fill_manual(COLOR_BLIND_PALETTE)
            )
    p.save(plotfile, height=2.7, width=(1.2 + 0.25 * len(names)),
            verbose=False)
    plt.close()


def plotReadsPerBC(names, readsperbcfiles, plotfile, 
        maxreads=10, maxcol=6):
    """Plots ``dms2_bcsubamp`` reads-per-barcode stats for set of samples.

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
    assert len(names) == len(readsperbcfiles) == len(set(names))
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

    # make name a category to preserve order
    df['name'] = df['name'].astype(
            pandas.api.types.CategoricalDtype(categories=names))

    ncol = min(maxcol, len(names))
    nrow = math.ceil(len(names) / float(ncol))
    p = (ggplot(df)
            + geom_col(aes(x='number of reads', y='number of barcodes'),
                position='stack')
            + scale_x_continuous(breaks=[1, maxreads // 2, maxreads],
                    labels=['$1$', '${0}$'.format(maxreads // 2), 
                    '$\geq {0}$'.format(maxreads)])
            + scale_y_continuous(labels=latexSciNot)
            + facet_wrap('~name', ncol=ncol)
            + theme(figure_size=(1.9 * (0.8 + ncol), 1.3 * (0.4 + nrow)))
            )
    p.save(plotfile, verbose=False)
    plt.close()


def plotDepth(names, countsfiles, plotfile, maxcol=4, charlist=CODONS):
    """Plot sequencing depth along primary sequence.

    Args:
        `names` (list or series)
            Names of samples for which we plot statistics.
        `countsfiles` (list or series)
            Files containing character counts at each site.
            Should have column named `site` and a column
            for each character in `charlist`.
        `plotfile` (str)
            Name of created PDF plot file containing count depth.
        `maxcol` (int)
            Number of columns in faceted plot.
        `charlist` (list)
            Characters contained in `countsfiles`. For instance,
            list of codons or amino acids.
    """
    assert len(names) == len(countsfiles) == len(set(names))
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'

    counts = pandas.concat(
            [pandas.read_csv(f)
                   .assign(name=name)
                   .assign(ncounts=lambda x: x[charlist].sum(axis=1))
                   .rename(columns={'ncounts':'number of counts'})
            for (name, f) in zip(names, countsfiles)], ignore_index=True)
    
    ncol = min(maxcol, len(names))
    nrow = math.ceil(len(names) / float(ncol))

    # make name a category to preserve order
    counts['name'] = counts['name'].astype(
            pandas.api.types.CategoricalDtype(categories=names))

    p = (ggplot(counts, aes(x='site', y='number of counts'))
            + geom_step(size=0.4)
            + scale_y_continuous(labels=latexSciNot, 
                    limits=(0, counts['number of counts'].max()))
            + scale_x_continuous(limits=(counts['site'].min(),
                    counts['site'].max()))
            + facet_wrap('~name', ncol=ncol) 
            + theme(figure_size=(2.25 * (0.6 + ncol), 1.3 * (0.3 + nrow)))
            )
    p.save(plotfile, verbose=False)
    plt.close()


def plotMutFreq(names, countsfiles, plotfile, maxcol=4):
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
    assert len(names) == len(countsfiles) == len(set(names))
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'

    counts = pandas.concat([dms_tools2.utils.annotateCodonCounts(f).assign(
            name=name) for (name, f) in zip(names, countsfiles)], 
            ignore_index=True).rename(columns={'mutfreq':'mutation frequency'})
    
    ncol = min(maxcol, len(names))
    nrow = math.ceil(len(names) / float(ncol))

    # make name a category to preserve order
    counts['name'] = counts['name'].astype(
            pandas.api.types.CategoricalDtype(categories=names))

    p = (ggplot(counts, aes(x='site', y='mutation frequency'))
            + geom_step(size=0.4)
            + scale_y_continuous(labels=latexSciNot, 
                    limits=(0, counts['mutation frequency'].max()))
            + scale_x_continuous(limits=(counts['site'].min(),
                    counts['site'].max()))
            + facet_wrap('~name', ncol=ncol) 
            + theme(figure_size=(2.25 * (0.6 + ncol), 1.3 * (0.3 + nrow)))
            )
    p.save(plotfile, verbose=False)
    plt.close()


def plotCumulMutCounts(names, countsfiles, plotfile, chartype,
        nmax=15, maxcol=4):
    """Plot fraction of mutations seen <= some number of times.

    For each set of counts in `countsfiles`, plot the fraction
    of mutations seen greater than or equal to some number of
    times. This is essentially a cumulative fraction plot.

    Args:
        `names` (list or series)
            Names of samples for which we plot statistics.
        `countsfiles` (list or series)
            ``*_codoncounts.csv`` files of type created by ``dms2_bcsubamp``.
        `plotfile` (str)
            Name of created PDF plot file.
        `chartype` (str)
            The type of character in `countsfiles`.

                - `codon` 

        `nmax` (int)
            Plot out to this number of mutation occurrences.
        `maxcol` (int)
            Number of columns in faceted plot.
    """
    assert len(names) == len(countsfiles) == len(set(names))
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'

    counts = pandas.concat([pandas.read_csv(f).assign(name=name) for
            (name, f) in zip(names, countsfiles)], ignore_index=True)

    if chartype != 'codon':
        raise ValueError("invalid chartype of {0}".format(chartype))

    codoncounts = pandas.concat([pandas.read_csv(f).assign(name=name)
            for (name, f) in zip(names, countsfiles)],
            ignore_index=True).assign(character='codons')
    assert set(CODONS) <= set(codoncounts.columns)
    codonmelt = codoncounts.melt(id_vars=['name', 'wildtype', 'character'], 
            value_vars=CODONS, value_name='counts',
            var_name='codon')
    codonmelt = codonmelt[codonmelt['codon'] != codonmelt['wildtype']]

    aacounts = pandas.concat([dms_tools2.utils.codonToAACounts(
            pandas.read_csv(f)).assign(name=name)
            for (name, f) in zip(names, countsfiles)],
            ignore_index=True).assign(character='amino acids')
    assert set(AAS_WITHSTOP) <= set(aacounts.columns)
    aamelt = aacounts.melt(id_vars=['name', 'character', 'wildtype'], 
            value_vars=AAS_WITHSTOP, value_name='counts',
            var_name='aa')
    aamelt = aamelt[aamelt['aa'] != aamelt['wildtype']]

    df = pandas.concat([codonmelt, aamelt], ignore_index=True)

    # make name a category to preserve order
    df['name'] = df['name'].astype(
            pandas.api.types.CategoricalDtype(categories=names))

    ncol = min(maxcol, len(names))
    nrow = math.ceil(len(names) / float(ncol))
    p = (ggplot(df, aes('counts', color='character', linestyle='character'))
            + stat_ecdf(geom='step', size=1)
            + coord_cartesian(xlim=(0, nmax))
            + facet_wrap('~name', ncol=ncol) 
            + theme(figure_size=(2.25 * (0.6 + ncol), 1.3 * (0.5 + nrow)),
                    legend_position='top', legend_direction='horizontal')
            + labs(color="") 
            + guides(color=guide_legend(title_position='left'))
            + ylab('fraction $\leq$ this many counts')
            + scale_color_manual(COLOR_BLIND_PALETTE)
            )
    p.save(plotfile, verbose=False)
    plt.close()


def plotCodonMutTypes(names, countsfiles, plotfile,
        classification='aachange', csvfile=None):
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

                `aachange` : stop, synonymous, nonsynonymous

                `n_ntchanges` : number of nucleotide changes per codon

                `singlentchanges` : nucleotide change in 1-nt mutations

        `csvfile` (str or `None`)
            `None` or name of CSV file to which numerical data are written.
    """
    assert len(names) == len(countsfiles) == len(set(names))
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
    elif classification == 'singlentchanges':
        muttypes = dict([(ntchange, ntchange) for ntchange in [
                '{0}to{1}'.format(nt1, nt2) for nt1 in NTS
                for nt2 in NTS if nt1 != nt2]])
    else:
        raise ValueError("Invalid classification {0}".format(classification))

    df = (counts[list(muttypes.values()) + ['ncounts', 'name']]
            .groupby('name', as_index=False)
            .sum(axis=1)
            .assign(ncounts=lambda x: x['ncounts'].astype('float'))
            )
    for (newcol, n) in muttypes.items():
        df[newcol] = (df[n] / df['ncounts']).fillna(0)

    if csvfile:
        df[['name'] + list(muttypes.keys())].to_csv(csvfile, index=False) 

    df = df.melt(id_vars='name', var_name='mutation type',
                value_vars=list(muttypes.keys()),
                value_name='per-codon frequency')

    p = (ggplot(df)
            + geom_col(aes(x='name', y='per-codon frequency', 
                fill='mutation type'), position='stack')
            + theme(axis_text_x=element_text(angle=90, vjust=1, hjust=0.5),
                    axis_title_x=element_blank())
            + scale_y_continuous(labels=latexSciNot)
            )
    if len(muttypes) <= len(COLOR_BLIND_PALETTE):
        p = p + scale_fill_manual(COLOR_BLIND_PALETTE)
    else:
        p = p + guides(fill=guide_legend(ncol=2))

    p.save(plotfile, height=2.7, width=(1.2 + 0.25 * len(names)),
            verbose=False)
    plt.close()


def plotCorrMatrix(names, infiles, plotfile, datatype,
        trim_unshared=True, title='', colors='black',
        contour=False, ncontours=10):
    """Plots correlations among replicates.

    Args:
        `names` (list or series)
            Names of samples for which we plot statistics.
        `infiles` (list or series)
            CSV files containing data. Format depends on `datatype`.
        `plotfile` (str)
            Name of created PDF plot file.
        `datatype` (str)
            Type of data for which we are plotting correlations:
                - `prefs`: in format returned by ``dms2_prefs``
                - `mutdiffsel`: mutdiffsel from ``dms2_diffsel``
                - `abs_diffsel`: sitediffsel from ``dms2_diffsel``
                - `positive_diffsel`: sitediffsel from ``dms2_diffsel``
                - `max_diffsel`: sitediffsel from ``dms2_diffsel``
                - `mutfracsurvive`: from ``dms2_fracsurvive``
        `trim_unshared` (bool)
            What if files in `infiles` don't have same sites / mutations?
            If `True`, trim unshared one and just analyze ones
            shared among all files. If `False`, raise an error.
        `title` (str)
            Title to place above plot.
        `colors` (str or list)
            Color(s) to color scatter points. If a string, should
            specify one color for all plots. Otherwise should be
            list of length `len(names) * (len(names) - 1) // 2`
            giving lists of colors for plots from top to bottom, 
            left to right.
        `contour` (bool)
            Show contour lines from KDE rather than points.
        `ncontours` (int)
            Number of contour lines if using `contour`.
    """
    assert len(names) == len(infiles) == len(set(names)) > 1
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'

    if datatype == 'prefs':
        # read prefs into dataframe, ensuring all have same characters
        prefs = [pandas.read_csv(f).assign(name=name) for (name, f) 
                in zip(names, infiles)]
        chars = set(prefs[0].columns)
        sites = set(prefs[0]['site'].values)
        for p in prefs:
            unsharedchars = chars.symmetric_difference(set(p.columns))
            if unsharedchars:
                raise ValueError("infiles don't have same characters: {0}"
                        .format(unsharedchars))
            unsharedsites = sites.symmetric_difference(set(p['site']))
            if trim_unshared:
                sites -= unsharedsites
            elif unsharedsites:
                raise ValueError("infiles don't have same sites: {0}".
                        format(unsharedsites))
        assert {'site', 'name'} < chars
        chars = list(chars - {'site', 'name'})
        # get measurements for each replicate in its own column
        df = (pandas.concat(prefs, ignore_index=True)
                    .query('site in @sites') # only keep shared sites
                    .melt(id_vars=['name', 'site'], var_name='char')
                    .pivot_table(index=['site', 'char'], columns='name')
                    )
        df.columns = df.columns.get_level_values(1)

    elif datatype in ['mutdiffsel', 'mutfracsurvive']:
        mut_df = [pandas.read_csv(f)
                        .assign(name=name)
                        .assign(mutname=lambda x: x.wildtype +
                                x.site.map(str) + x.mutation)
                        .sort_values('mutname')
                        [['name', 'mutname', datatype]]
                    for (name, f) in zip(names, infiles)]
        muts = set(mut_df[0]['mutname'].values)
        for m in mut_df:
            unsharedmuts = muts.symmetric_difference(set(m['mutname']))
            if trim_unshared:
                muts -= unsharedmuts
            elif unsharedmuts:
                raise ValueError("infiles don't have same muts: {0}".
                        format(unsharedmuts))
        df = (pandas.concat(mut_df, ignore_index=True)
                    .query('mutname in @muts') # only keep shared muts
                    .pivot_table(index='mutname', columns='name')
                    .dropna()
                    )
        df.columns = df.columns.get_level_values(1)

    elif datatype in ['abs_diffsel', 'positive_diffsel', 'max_diffsel',
            'avgfracsurvive', 'maxfracsurvive']:
        site_df = [pandas.read_csv(f)
                         .assign(name=name)
                         .sort_values('site')
                         [['name', 'site', datatype]]
                    for (name, f) in zip(names, infiles)]
        sites = set(site_df[0]['site'].values)
        for s in site_df:
            unsharedsites = sites.symmetric_difference(set(s['site']))
            if trim_unshared:
                sites -= unsharedsites
            elif unsharedsites:
                raise ValueError("infiles don't have same sites: {0}".
                        format(unsharedsites))
        df = (pandas.concat(site_df, ignore_index=True)
                    .query('site in @sites') # only keep shared sites
                    .pivot_table(index='site', columns='name')
                    .dropna()
                    )
        df.columns = df.columns.get_level_values(1)

    else:
        raise ValueError("Invalid datatype {0}".format(datatype))

    ncolors = len(names) * (len(names) - 1) // 2
    if isinstance(colors, str):
        colors = [colors] * ncolors
    else:
        assert len(colors) == ncolors, "not {0} colors".format(ncolors)

    def corrfunc(x, y, contour, **kws):
        r, _ = scipy.stats.pearsonr(x, y)
        ax = plt.gca()
        ax.annotate('R = {0:.2f}'.format(r), xy=(0.05, 0.9), 
                xycoords=ax.transAxes, fontsize=19, 
                fontstyle='oblique')
        color = colors.pop(0)
        if contour:
            seaborn.kdeplot(x, y, shade=True, n_levels=ncontours)
        else:
            plt.scatter(x, y, s=22, alpha=0.35, color=color,
                    marker='o', edgecolor='none', rasterized=True)

    # map lower / upper / diagonal as here:
    # https://stackoverflow.com/a/30942817 for plot
    p = seaborn.PairGrid(df)
    p.map_lower(corrfunc, colors=colors, contour=contour)
    if datatype == 'prefs':
        p.set(  xlim=(0, 1), 
                ylim=(0, 1), 
                xticks=[0, 0.5, 1],
                yticks=[0, 0.5, 1],
                xticklabels=['0', '0.5', '1'],
                yticklabels=['0', '0.5', '1'],
                )
        # hide upper, diag: https://stackoverflow.com/a/34091733
        for (i, j) in zip(*numpy.triu_indices_from(p.axes, 1)):
            p.axes[i, j].set_visible(False)
        for (i, j) in zip(*numpy.diag_indices_from(p.axes)):
            p.axes[i, j].set_visible(False)

    elif datatype in ['mutdiffsel', 'abs_diffsel', 'positive_diffsel',
            'max_diffsel', 'mutfracsurvive', 'avgfracsurvive',
            'maxfracsurvive']:

        p.map_diag(seaborn.distplot, color='black', kde=True, hist=False)

        (lowlim, highlim) = (df.values.min(), df.values.max())
        highlim += (highlim - lowlim) * 0.05
        lowlim -= (highlim - lowlim) * 0.05
        p.set(xlim=(lowlim, highlim), ylim=(lowlim, highlim))

        # hide upper
        for (i, j) in zip(*numpy.triu_indices_from(p.axes, 1)):
            p.axes[i, j].set_visible(False)

    else:
        raise ValueError("invalid datatype")
        p.map_upper(seaborn.kdeplot, n_levels=20, cmap='Blues_d')

    if title:
        # following here: https://stackoverflow.com/a/29814281
        p.fig.suptitle(title)

    p.savefig(plotfile)
    plt.close()


def plotSiteDiffSel(names, diffselfiles, plotfile, 
        diffseltype, maxcol=2, white_bg=False):
    """Plot site diffsel or fracsurvive along sequence.

    Despite the function name, this function can be used to
    plot either differential selection or fraction surviving.

    Args:
        `names` (list or series)
            Names of samples for which we plot statistics.
        `diffselfiles` (list or series)
            ``*sitediffsel.csv`` files from ``dms2_diffsel`` or
            ``*sitefracsurvive.csv`` files from ``dms2_fracsurvive``.
        `plotfile` (str)
            Name of created PDF plot file.
        `diffseltype` (str)
            Type of diffsel or fracsurvive to plot:
                - `positive`: positive sitediffsel
                - `total`: positive and negative sitediffsel
                - `max`: maximum mutdiffsel
                - `minmax`: minimum and maximum mutdiffsel
                - `avgfracsurvive`: total site fracsurvive
                - `maxfracsurvive`: max mutfracsurvive at site
        `maxcol` (int)
            Number of columns in faceted plot.
        `white_bg` (bool)
            Plots will have a white background with limited other formatting.

    """
    assert len(names) == len(diffselfiles) == len(set(names)) > 0
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'

    diffsels = [pandas.read_csv(f).assign(name=name) for (name, f) 
            in zip(names, diffselfiles)]
    assert all([set(diffsels[0]['site']) == set(df['site']) for df in 
            diffsels]), "diffselfiles not all for same sites"
    diffsel = pandas.concat(diffsels, ignore_index=True)

    ylabel = 'differential selection'
    if diffseltype == 'positive':
        rename = {'positive_diffsel':'above'}
    elif diffseltype == 'total':
        rename = {'positive_diffsel':'above',
                  'negative_diffsel':'below'}
    elif diffseltype == 'max':
        rename = {'max_diffsel':'above'}
    elif diffseltype == 'minmax':
        rename = {'max_diffsel':'above',
                  'min_diffsel':'below'}
    elif diffseltype in ['avgfracsurvive', 'maxfracsurvive']:
        ylabel = 'fraction surviving'
        rename = {diffseltype:'above'}
    else:
        raise ValueError("invalid diffseltype {0}".format(diffseltype))
    diffsel = (diffsel.rename(columns=rename)
                      .melt(id_vars=['site', 'name'], 
                            value_vars=list(rename.values()),
                            value_name='diffsel',
                            var_name='direction')
                      )


    # natural sort by site: https://stackoverflow.com/a/29582718
    diffsel = diffsel.reindex(index=natsort.order_by_index(
            diffsel.index, natsort.index_natsorted(diffsel.site,
            signed=True)))
    # now some manipulations to make site str while siteindex is int
    diffsel['site'] = diffsel['site'].apply(str)
    diffsel['siteindex'] = pandas.Categorical(diffsel['site'],
            diffsel['site'].unique()).codes
    
    ncol = min(maxcol, len(names))
    nrow = math.ceil(len(names) / float(ncol))

    # make name a category to preserve order
    diffsel['name'] = diffsel['name'].astype(
            pandas.api.types.CategoricalDtype(categories=names))

    (xbreaks, xlabels) = breaksAndLabels(diffsel['siteindex'].unique(), 
            diffsel['site'].unique(), n=6)
    if white_bg:
        p = (ggplot(diffsel, aes(x='siteindex', y='diffsel',
                    color='direction', fill='direction'))
             + geom_step(size=0.3)
             + xlab('site')
             + ylab(ylabel)
             + scale_x_continuous(breaks=xbreaks, labels=xlabels)
             + scale_color_manual(COLOR_BLIND_PALETTE)
             + scale_fill_manual(COLOR_BLIND_PALETTE)
             + guides(color=False)
             + theme(panel_background=element_rect(fill='white'),
                     axis_line_x=element_line(color='black'),
                     axis_line_y=element_line(color='black'),
                     panel_grid=element_blank(),
                     panel_border=element_blank(),
                     strip_background=element_blank()
                     )
             )
    else:
        p = (ggplot(diffsel, aes(x='siteindex', y='diffsel', color='direction'))
             + geom_step(size=0.4)
             + xlab('site')
             + ylab(ylabel)
             + scale_x_continuous(breaks=xbreaks, labels=xlabels)
             + scale_color_manual(COLOR_BLIND_PALETTE)
             + guides(color=False)
             )
    if not ((len(names) == 1) and ((not names[0]) or names[0].isspace())):
        p = p + facet_wrap('~name', ncol=ncol)
    p = p + theme(figure_size=(4.6 * (0.3 + ncol), 1.9 * (0.2 + nrow)))
    p.save(plotfile, verbose=False)
    plt.close()


def plotFacetedNeutCurves(
        neutdata,
        plotfile,
        xlabel,
        ylabel,
        maxcol=3):
    """Faceted neutralization curves with points and fit line.

    Args:
        `neutdata` (pandas DataFrame)
            Should have the following columns:
            `concentration`, `sample`, `fit`, `points`.
            The plot is faceted on `sample`. The line
            smoothly connects all points in column
            `fit`, and points are drawn anywhere
            that `points` is not `NaN`.
        `plotfile` (str)
            Name of created plot.
        `xlabel` (str)
            x-axis label
        `ylabel` (str)
            y-axis label
        `maxcol` (int)
            Number of columns in facets.
    """
    cols = {'concentration', 'sample', 'fit', 'points'}
    assert set(neutdata.columns) >= cols, ("missing cols:\n"
            "required: {0}\nactual: {1}".format(cols, neutdata.columns))

    # make sample a category to preserve order
    neutdata = neutdata.copy()
    samples = neutdata['sample'].unique()
    neutdata['sample'] = neutdata['sample'].astype(
            pandas.api.types.CategoricalDtype(categories=samples))

    ncol = min(maxcol, len(samples))
    nrow = math.ceil(len(samples) / float(ncol))

    ymin = min(neutdata['fit'].min(), neutdata['points'].min(), 0)
    ymax = max(neutdata['fit'].max(), neutdata['points'].max(), 0)

    p = (ggplot(neutdata) +
            geom_point(aes(x='concentration', y='points')) +
            geom_line(aes(x='concentration', y='fit')) +
            scale_x_log10(labels=latexSciNot) +
            scale_y_continuous(limits=(ymin, ymax)) +
            xlab(xlabel) + 
            ylab(ylabel) +
            facet_wrap('~sample', ncol=ncol) 
            + theme(figure_size=(2.4 * (0.25 + ncol),
                                 1.45 * (0.25 + nrow)))
            )
    p.save(plotfile, verbose=False)
    plt.close()


def findSigSel(df, valcol, plotfile, fdr=0.05, title=None):
    """Finds "significant" selection at sites / mutations.

    Designed for the case where most sites / mutations are not
    under selection, but a few may be. It tries to find those
    few that are under selection. 

    It does not use a mechanistic statistical model, but rather uses
    `robust regression <http://scipy-cookbook.readthedocs.io/items/robust_regression.html>`_
    (soft L1 loss) to fit a gamma distribution. The rationale for
    a gamma distribution is that it is negative binomial's continuous
    `analog <http://www.nehalemlabs.net/prototype/blog/2013/12/01/gamma-distribution-approximation-to-the-negative-binomial-distribution/>`_.
    It then identifies sites that clearly have **larger** values than
    expected under this distribution. It currently does not
    identify sites with **smaller** (or more negative) than expected
    values.

    Args:
        `df` (pandas DataFrame)
            Contains data to analyze
        `valcol` (string)
            Column in `df` with values (e.g., `fracsurvive`)
        `plotfile` (string)
            Name of file to which we plot fit.
        `fdr` (float)
            Find sites that are significant at this `fdr`
            given fitted distribution.
        `title` (string or `None`)
            Title for plot.

    Returns:
        Creates the plot in `plotfile`. Also returns
        the 3-tuple `(df_sigsel, cutoff, gamma_fit)` where:
        
            - `df_sigsel` is copy of `df` with new columns
              `P`, `Q`, and `sig`. These give P value, Q
              value, and whether site meets `fdr` cutoff
              for significance.

            - `cutoff` is the maximum value that is **not**
              called significant. Because FDR is a property
              of a distribution, this value cannot be 
              interpreted as meaning a new data point would
              be called based on this cutoff, as the cutoff
              would change. But `cutoff` is useful for 
              plotting to see significant / non-significant.

            - `gamma_params` is a `numpy.ndarray` that of 
              length 3 that gives the shape, scale, and location
              parameter of the fit gamma distribution.

    An example: First, simulate points from a gamma distribution:

    >>> shape_sim = 1.5
    >>> scale_sim = 0.005
    >>> loc_sim = 0.0
    >>> gamma_sim = scipy.stats.gamma(shape_sim, scale=scale_sim,
    ...         loc=loc_sim)
    >>> nsites = 1000
    >>> scipy.random.seed(0)
    >>> df = pandas.DataFrame.from_dict({
    ...         'site':[r for r in range(nsites)],
    ...         'fracsurvive':gamma_sim.rvs(nsites)})

    Now make two sites have "significantly" higher values:

    >>> sigsites = [100, 200]
    >>> df.loc[sigsites, 'fracsurvive'] = 0.08

    Now plot and find the significant sites:

    >>> plotfile = '_findSigSel.png'
    >>> (df_sigsel, cutoff, gamma_params) = findSigSel(
    ...         df, 'fracsurvive', plotfile, title='example')

    Here is the resulting plot:

    .. image:: _static/_findSigSel.png
       :width: 4in
       :align: center

    Make sure the fitted params are close to the ones used to
    simulate the data:

    >>> numpy.allclose(shape_sim, gamma_params[0], rtol=0.1, atol=1e-3)
    True
    >>> numpy.allclose(scale_sim, gamma_params[1], rtol=0.1, atol=1e-3)
    True
    >>> numpy.allclose(loc_sim, gamma_params[2], rtol=0.1, atol=1e-3)
    True

    Check that we find the correct significant sites:

    >>> set(sigsites) == set(df_sigsel.query('sig').site)
    True

    Make sure that sites above cutoff are significant:

    >>> df_sigsel.query('sig').equals(df_sigsel.query('fracsurvive > @cutoff'))
    True
    """
    assert valcol in df.columns, "no `valcol` {0}".format(valcol)

    newcols = {'P', 'Q', 'sig'}
    assert not (newcols & set(df.columns)), \
            "`df` already has {0}".format(newcols)

    def _f(x, bins, heights):
        """Gamma distribution least squares fitting function.

        Zero when distribution perfectly fits histogram.
        `x` is `(shape, scale, loc)`.
        """
        return (scipy.stats.gamma.pdf(bins, x[0], scale=x[1],
                loc=x[2]) - heights)

    # We fit curves to histogram. First we need to get bins.
    try:
        # try with Freedman Diaconis Estimator
        binedges = numpy.histogram(df[valcol], bins='fd')[1]
    except ValueError:
        # fd will fail of lots of identical points
        binedges = numpy.histogram(df[valcol], bins='doane')[1]

    # get bin centers
    bins = (binedges[ : -1] + binedges[1 : ]) / 2

    # plot the histogram
    plt.figure(figsize=(5.5, 4))
    (heights, binedges, patches) = plt.hist(df[valcol],
            bins=binedges, normed=True, histtype='stepfilled',
            color=COLOR_BLIND_PALETTE[2])

    # initial guess gives correct mean and variance for
    # gamma distribution with loc of 0
    scale = df[valcol].var() / df[valcol].mean()
    shape = df[valcol].mean() / scale
    x0 = numpy.array([shape, scale, 0.0])

    # fit using soft L1 loss for robust regression
    # http://scipy-cookbook.readthedocs.io/items/robust_regression.html
    fit = scipy.optimize.least_squares(_f, x0, args=(bins, heights),
            loss='soft_l1')
    gamma_params = fit.x
    gamma_fit = scipy.stats.gamma(fit.x[0], scale=fit.x[1],
            loc=fit.x[2])

    # add fit gamma distribution to plot
    nfitbins = 500
    if nfitbins > len(bins):
        fitbins = numpy.linspace(bins[0], bins[-1], nfitbins)
    else:
        fitbins = bins
    plt.plot(fitbins, gamma_fit.pdf(fitbins),
            color=COLOR_BLIND_PALETTE[1])

    # compute P and Q values
    df_sigsel = (df.assign(P=lambda x: gamma_fit.sf(x[valcol]))
                   .assign(Q=lambda x: multipletests(x.P, fdr, 'fdr_bh')[1])
                   .assign(sig=lambda x: multipletests(x.P, fdr, 'fdr_bh')[0])
                   )

    # compute cutoff 
    cutoff = df_sigsel.query('not sig')[valcol].max()

    # plot cutoff
    # find first bin boundary greater than cutoff
    if (binedges > cutoff).any():
        bincutoff = binedges[binedges > cutoff][0]
    else:
        bincutoff = binedges[-1]
    # now annotate plot
    plt.axvline(bincutoff, color=COLOR_BLIND_PALETTE[3], ls='--', lw=0.75)
    text_y = 0.95 * plt.ylim()[1]
    (xmin, xmax) = plt.xlim()
    if (bincutoff - xmin) < 0.75 * (xmax - xmin):
        text_x = bincutoff + 0.01 * (xmax - xmin)
        ha = 'left'
    else:
        text_x = bincutoff - 0.01 * (xmax - xmin)
        ha = 'right'
    if len(df_sigsel.query('sig')):
        text = '{0} values\nsignificant\n($>${1})'.format(
                len(df_sigsel.query('sig')), latexSciNot([cutoff])[0])
    else:
        text = 'no values\nsignificant'
    plt.text(text_x, text_y, text,
            horizontalalignment=ha, verticalalignment='top',
            color=COLOR_BLIND_PALETTE[3], size='small')

    # put labels on plot
    plt.xlabel(valcol.replace('_', ' '))
    plt.ylabel('density')
    if title:
        plt.title(title.replace('_', ' '))

    # save plot
    plt.tight_layout()
    plt.savefig(plotfile)
    plt.close()

    return (df_sigsel, cutoff, gamma_params)


def plotColCorrs(df, plotfile, cols, *, lower_filter=None,
        title=None, shrink_threshold=25):
    """Plots correlation among columns in pandas Data Frame.

    Plots distribution of each variable and pairwise correlations.

    Args:
        `df` (pandas DataFrame)
            Data frame with data to plot.
        `plotfile` (str)
            Name of created plot.
        `cols` (list)
            List of columns in `df` to plot.
        `lower_filter` (`None` or str)
            Can be any string that can passed to the `query` function
            of `df`. In this case, on the lower diagonal only plot
            data for which this query is `True`.
        `title` (`None` or str)
            Title of plot.
        `shrink_threshold` (float)
            See argument of same name to :meth:`hist_bins_intsafe`.
    """
    if not set(cols).issubset(set(df.columns)):
        raise ValueError("`cols` specifies columns not in `df`")

    if lower_filter is not None:
        filter_indices = df.query(lower_filter).index
        color_all = COLOR_BLIND_PALETTE_GRAY[0]
        color_filter = COLOR_BLIND_PALETTE_GRAY[1]
    else:
        color_all = COLOR_BLIND_PALETTE[0]

    def hist1d(x, color, **kwargs):
        """1D histogram for diagonal elements."""
        bins=dms_tools2.plot.hist_bins_intsafe(x,
                shrink_threshold=shrink_threshold)
        plt.hist(x, color=color_all, bins=bins, **kwargs)
        if lower_filter:
            plt.hist(x.ix[filter_indices], color=color_filter,
                     bins=bins, **kwargs)

    def hist2d(x, y, color, filterdata, **kwargs):
        """2D histogram for off-diagonal elements."""
        bins = [dms_tools2.plot.hist_bins_intsafe(a,
                shrink_threshold=shrink_threshold) for a in [x, y]]
        if filterdata:
            color = color_filter
            x = x.ix[filter_indices]
            y = y.ix[filter_indices]
        else:
            color = color_all
        cmap = dms_tools2.plot.from_white_cmap(color)
        plt.hist2d(x, y, bins=bins, cmap=cmap, **kwargs)

    g = (dms_tools2.plot.AugmentedPairGrid(df, vars=cols,
            diag_sharey=False, height=3)
         .map_diag(hist1d)
         .map_upper(hist2d, filterdata=False)
         .map_lower(hist2d, filterdata=(lower_filter is not None))
         .ax_lims_clip_outliers()
         )

    if lower_filter is not None:
        label_order = ['all\n({0})'.format(
                            dms_tools2.plot.latexSciNot(len(df))),
                       '{0}\n({1})'.format(lower_filter,
                            dms_tools2.plot.latexSciNot(
                            len(df.query(lower_filter)))),
                       ]
        label_data = {lab:plt.Line2D([0], [0], color=c, lw=10,
                                     solid_capstyle='butt')
                      for (lab, c) in zip(label_order,
                          [color_all, color_filter])
                      }
        g.add_legend(label_data, label_order=label_order,
                     labelspacing=2, handlelength=1.5)

    if title is not None:
        g.fig.suptitle(title, va='bottom')

    g.savefig(plotfile)
    plt.close()


def plotRarefactionCurves(df, rarefy_col, plotfile,
        *, facet_col=None, nrow=1,
        xlabel='reads', ylabel=None, facet_scales='free'):
    """Plots rarefaction curves.

    The rarefaction curves are calculated analytically using
    :py:mod:`dms_tools2.utils.rarefactionCurve`.

    Args:
        `df` (pandas DataFrame)
            Data frame containing data. In tidy form if
            faceting.
        `rarefy_col` (str)
            Name of column in `df` that contains the variable
            that we rarify. For instance, these might be strings
            giving barcodes.
        `plotfile` (str)
            Name of created plot.
        `facet_col` (str or `None`)
            If not `None`, should be name of a column in `df`
            that contains a variable we facet in the plot.
        `nrow` (int)
            If faceting, the number of rows.
        `xlabel` (str)
            X-axis label.
        `ylabel` (str or `None`)
            Y-axis label. If `None`, defaults to value of `rarefy_col`.
        `facet_scales` (str`)
            Scales for faceting. Can be "free", "free_x", "free_y", or
            "fixed"

    Here is an example. First, we simulate two sets of barcodes.
    For ease of fast simulation, the barcodes are just numbers here.
    One samples a large set, the other a set a quarter that size with
    half as many reads:

    >>> nbc = 40000
    >>> bclen = 10
    >>> nreads = 200000
    >>> barcodes = list(range(nbc))
    >>> numpy.random.seed(1)
    >>> large_set = numpy.random.choice(barcodes, size=nreads)
    >>> small_set = numpy.random.choice(barcodes[ : nbc // 4], size=nreads // 2)

    Now we put these in a tidy data frame where one column is named
    "barcodes" and the other is named "sample":

    >>> df = pandas.DataFrame({
    ...         "barcodes":list(small_set) + list(large_set),
    ...         "sample":['small_set'] * len(small_set) +
    ...                  ['large_set'] * len(large_set)})

    Finally, plot the rarefaction curves:

    >>> plotfile = '_plotRarefactionCurves.png'
    >>> plotRarefactionCurves(df, 'barcodes', plotfile, facet_col='sample')

    Here is the resulting plot:

    .. image:: _static/_plotRarefactionCurves.png
       :width: 6in
       :align: center
    """
    if rarefy_col not in df.columns:
        raise ValueError("`df` does not have `rarefy_col` {0}"
                         .format(rarefy_col))
    if (facet_col is not None) and (facet_col not in df.columns):
        raise ValueError("`df` does not have `facet_col` {0}"
                         .format(facet_col))
    if ylabel is None:
        ylabel = rarefy_col

    # get iterator over groups or dummy iterator
    categories = None
    if facet_col is not None:
        if df[facet_col].dtype.name == 'category':
            categories = df[facet_col].cat.categories
        df_iterator = df.groupby(facet_col)[rarefy_col]
        nfacets = len(df[facet_col].unique())
    else:
        df_iterator = ('dummy', df[rarefy_col])
        nfacets = 1
    d = collections.defaultdict(list)

    for name, group in df_iterator:
        xs, ys = dms_tools2.utils.rarefactionCurve(group)
        assert len(xs) == len(ys)
        d[ylabel] += ys
        d[xlabel] += xs
        d['_facet_var'] += [name] * len(xs)
    rarefied = pandas.DataFrame(d)
    if categories is not None:
        rarefied['_facet_var'] = pandas.Categorical(
                rarefied['_facet_var'], categories)

    ident = lambda x: x.astype('int') if all(x.astype('int') == x) else x
    if rarefied[xlabel].max() >= 1e4: 
        xlabeler = latexSciNot
    else:
        xlabeler = ident
    if rarefied[ylabel].max() >= 1e4:
        ylabeler = latexSciNot
    else:
        ylabeler = ident

    p = (ggplot(rarefied, aes(xlabel, ylabel)) +
            geom_line() +
            xlab(xlabel) +
            ylab(ylabel) +
            scale_x_continuous(labels=xlabeler,
                    breaks=lambda x: matplotlib.ticker.MaxNLocator(3)
                                     .tick_values(min(x), max(x))) +
            scale_y_continuous(labels=ylabeler)
            )

    x_panel_spacing = {"free":0.75, "free_x":0.1,
                       "fixed":0.1, "free_y":0.75}[facet_scales]
    y_panel_spacing = {"free":0.4, "free_x":0.4,
                       "fixed":0.1, "free_y":0.1}[facet_scales]
    if facet_col is not None:
        p = p + facet_wrap('~ _facet_var', nrow=nrow, scales=facet_scales)
        p = p + theme(panel_spacing_x=x_panel_spacing,
                   panel_spacing_y=y_panel_spacing)

    ncol = math.ceil(nfacets / nrow)
    p.save(plotfile,
            height=0.5 + 1.75 * nrow + y_panel_spacing * (nrow - 1),
            width=(1.25 + 2 * ncol + x_panel_spacing * (ncol - 1)),
            verbose=False,
            limitsize=False)
    plt.close()


def hist_bins_intsafe(x, method='fd', shrink_threshold=None, maxbins=100):
    """Histogram bins that work for integer data.

    You can auto-choose bins using `numpy.histogram`.
    However, if the data are integer, these bins
    may be non-integer and so some bins will capture
    more integers. This function fixes that.

    Args:
        `x` (numpy array)
            The data to bin.
        `method` (str)
            The binning method. Can be anything
            acceptable to `numpy.histogram` as
            a `bins` argument.
        `shrink_threshold` (`None` or int)
            If set to a value other than `None`,
            apply a heuristic threshold to slow
            the growth in number of bins if they
            exceed this number.
        `maxbins` (int)
            Maximum number of bins.

    Returns:
        The bin edges as returned in the second
        element of `numpy.histogram`, but adjusted
        to be of integer width if the data are all
        integers.

    Just like `numpy.histogram` for non-int data:

    >>> numpy.random.seed(1)
    >>> x = 100 * numpy.random.random(500)
    >>> bin_edges = numpy.histogram(x, bins='fd')[1]
    >>> bin_edges_intsafe = hist_bins_intsafe(x)[ : len(bin_edges)]
    >>> numpy.allclose(bin_edges, bin_edges_intsafe)
    True
    >>> numpy.allclose(bin_edges, bin_edges.astype('int'))
    False

    But gives integer bins for int data:
    >>> x = x.astype('int')
    >>> bin_edges = numpy.histogram(x, bins='fd')[1]
    >>> bin_edges_intsafe = hist_bins_intsafe(x)
    >>> numpy.allclose(bin_edges, bin_edges_intsafe)
    False
    >>> numpy.allclose(bin_edges, bin_edges.astype('int'))
    False
    >>> numpy.allclose(bin_edges_intsafe, bin_edges_intsafe.astype('int'))
    True
    """
    bin_edges = numpy.histogram(x, bins='fd')[1]
    if len(bin_edges) > maxbins:
        bin_edges = numpy.histogram(x, min(maxbins, len(bin_edges) - 1))[1]
    if shrink_threshold is None:
        corr = 1
    else:
        assert shrink_threshold > 1
        corr = max(1,
                math.sqrt(len(bin_edges) / float(shrink_threshold)))

    binwidth = (bin_edges[1] - bin_edges[0]) * corr
    if (x.astype('int') == x).all():
        binwidth = math.ceil(binwidth)
    bin_edges = numpy.arange(x.min(), x.max() + binwidth, binwidth)
    return bin_edges


def from_white_cmap(color):
    """Get matplotlib color map from white to `color`."""
    light = seaborn.set_hls_values(color, l=1)
    return seaborn.blend_palette([light, color], None, True)


class AugmentedPairGrid(seaborn.PairGrid):
    """Augmented version of `seaborn.PairGrid`."""

    def ax_lims_clip_outliers(self, frac_clip=0.001, extend=0.03):
        """Sets axis limits to clip outliers in data.

        Useful if there are a few data points far outside the range
        of most of the data.

        Args:
            `frac_clip` (float)
                Set upper and lower limits so that this fraction of
                data is outside limits at both ends. Done **before**
                adding `extend`.
            `extend` (float)
                Extend the limits determined by `frac_clip` by this
                fraction of the data range.
        """
        assert 0 <= frac_clip < 0.5

        def _get_lims(s):
            """Gets limits for data in pandas.Series `s`."""
            s_range = s.max() - s.min()
            if s_range == 0:
                s_extend = max(extend, extend * s.max())
            else:
                s_extend = s_range * extend
            s_min = s.quantile(frac_clip) - s_extend
            s_max = s.quantile(1 - frac_clip) + s_extend
            return (s_min, s_max)

        xlims = [_get_lims(self.data[x]) for x in self.x_vars]
        ylims = [_get_lims(self.data[y]) for y in self.y_vars]

        for icol, xlim in enumerate(xlims):
            for irow, ylim in enumerate(ylims):
                self.axes[irow, icol].set_xlim(*xlim)
                if icol != irow:
                    self.axes[irow, icol].set_ylim(*ylim)

        return self


if __name__ == '__main__':
    import doctest
    doctest.testmod()
