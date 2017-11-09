"""
===================
rplot
===================

Plotting that uses ``R``.

The `dms_tools2` is written in Python, but there are
useful plotting features that are only available in ``R``,
such `ggseqlogo <https://omarwagih.github.io/ggseqlogo/>`_.

This module uses ``R`` to make plots using these features. It
requires `rpy2 <https://rpy2.readthedocs.io>`_ to be installed.
Installation of `rpy2 <https://rpy2.readthedocs.io>`_ is 
**not** automatic when you install `dms_tools2` unless
you have ``R`` and `rpy2 <https://rpy2.readthedocs.io>`_
installed, or install `dms_tools2` using::

    pip install dms_tools[rplot] --user

"""

import os
import io

import pandas

import phydmslib.weblogo

#: default colors for amino acid letters, by functional group
AA_COLORS_FG = phydmslib.weblogo.FunctionalGroupColorMapping()[1]

# import rpy2
try:
    import rpy2
except ImportError:
    raise ImportError("You must install `rpy2` to use `rplot`")
import rpy2.rinterface as rinterface
rinterface.initr()
from rpy2.robjects import pandas2ri, r
pandas2ri.activate()
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector, ListVector
from rpy2.rlike.container import TaggedList
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

# install necessary packages
utils = importr('utils')
utils.chooseCRANmirror(ind=1)
_packages = ['ggplot2', 'ggseqlogo']
utils.install_packages(StrVector(_packages))
for _package in _packages:
    importr(_package)


def versionInfo():
    """Returns string giving `rpy2` and ``R`` version info."""
    return "Using `rpy2` {0} with ``R`` {1}".format(
            rpy2.__version__, rinterface.R_VERSION_BUILD[1])


def siteSubsetGGSeqLogo(logodata, letters, plotfile, width, height,
        letter_colors=AA_COLORS_FG):
    """Creates one-row logo plot with subset of sites.

    Designed to show logo plot for a subset of sites. This
    is useful when you have data for many sites, but only
    want to look at a few of them. 

    Args:
        `logodata` (pandas DataFrame)
            Contains data to plot. Should have the columns
            `site`, `show`, and a column giving the height
            height of each letter in `letters`. Only sites
            where `show` is `True` are shown. Sites are 
            shown in the order they occur in this dataframe,
            with spaces every time there is an interspersed
            site with `show` being `False`. 
        `letters` (list)
            Letters for which we plot heights.
        `plotfile` (str)
            Name of created plot.
        `width` (float)
            Width of plot in inches.
        `height` (float)
            Height of plot in inches.
        `letter_colors` (dict)
            Values give color for every letter in `letters`.

    Here is an example that creates a plot for a subset of
    sites for two letters:

    >>> logodata = pandas.read_csv(io.StringIO(
    ...     '''site show    A    C
    ...        A101 True  0.8  0.2
    ...        N102 True  0.7  0.3
    ...        K103 False 0.1  0.9
    ...        L104 True  0.8  0.2
    ...        S105 True  0.5  0.5
    ...        T106 False 0.2  0.8
    ...        L107 True  0.7  0.3'''),
    ...     delim_whitespace=True, index_col=False)
    >>> plotfile = '_siteSubsetGGSeqLogo_test_plot.png'
    >>> siteSubsetGGSeqLogo(logodata,
    ...         letters=['A', 'C'],
    ...         plotfile=plotfile,
    ...         width=3, height=2.5
    ...         )
    >>> os.path.isfile(plotfile)
    True

    """
    if os.path.isfile(plotfile):
        os.remove(plotfile)

    assert set(letters) <= set(letter_colors.keys()), \
            "`letter_colors` not defined for all letters"

    expectcol = ['site', 'show'] + letters
    assert set(logodata.columns) <= set(expectcol), \
            "`logodata` needs these column: {0}".format(expectcol)

    # group consecutive rows into facets
    logodata = logodata[expectcol]
    logodata['facet'] = (
            (logodata['show'] != logodata['show'].shift(1))
            .astype(int).cumsum()
            )
    logodata = logodata.query('show')

    # generate list of matrices to facet
    matrices = []
    facets = logodata['facet'].unique()
    for f in facets:
        facetdata = (logodata.query('facet == @f')
                     .set_index('site')
                     [letters]
                     )
        m = r.matrix(
                facetdata.values.ravel(),
                ncol=len(facetdata),
                dimnames=[letters, facetdata.index.tolist()]
                )
        matrices.append(m)
    matrices = ListVector(TaggedList(matrices,
            tags=facets.astype('str')))

    print(matrices)



def facetedGGSeqLogo(logodata, letters, plotfile, width, height,
        ncol=None, letter_colors=AA_COLORS_FG):
    """Creates faceted logo plot.

    Designed to show several measurements on the same site
    site-by-side, potentially for many sites.

    Makes panel of logo plots faceted on `logodata['facetlabel']`,
    where letter stacks are labeled by `logodata['stacklabel']`
    and show the letters at the indicated heights.

    Args:
        `logodata` (pandas DataFrame)
            Contains data to plot. Should have the columns
            `facetlabel`, `stacklabel`, and a column giving the
            height of each letter in `letters`.
        `letters` (list)
            Letters for which we plot heights.
        `plotfile` (str)
            Name of created plot.
        `width` (float)
            Width of plot in inches.
        `height` (float)
            Height of plot in inches.
        `ncol` (int or `None`)
            Number of columns in faceted plot. If `None`, use
            as many as needed to plot everything in one row.
        `letter_colors` (dict)
            Values give color for every letter in `letters`.

    Here is an example that creates two facets each with
    two stacks for the letters `A` and `C`:

    >>> logodata = pandas.read_csv(io.StringIO(
    ...     '''facetlabel  stacklabel   A   C
    ...            site-1       BF520 0.8 0.2
    ...            site-1       BG505 0.9 0.1
    ...            site-2       BF520 0.4 0.6
    ...            site-2       BG505 0.5 0.5'''),
    ...     delim_whitespace=True, index_col=False)
    >>> plotfile = '_facetedGGSeqLogo_test_plot.png'
    >>> facetedGGSeqLogo(logodata,
    ...         letters=['A', 'C'],
    ...         plotfile=plotfile,
    ...         width=3, height=2.5
    ...         )
    >>> os.path.isfile(plotfile)
    True

    Here is the plot created by the code block above:

    .. image:: _static/_facetedGGSeqLogo_test_plot.png
       :width: 40%
       :align: center

    """
    if os.path.isfile(plotfile):
        os.remove(plotfile)

    assert set(letters) <= set(letter_colors.keys()), \
            "`letter_colors` not defined for all letters"

    # get and order data columns
    df_cols = ['facetlabel', 'stacklabel'] + letters
    assert set(logodata.columns) <= set(df_cols), "df lacks required columns"
    logodata = logodata[df_cols] 

    facets = logodata['facetlabel'].unique()
    stacks = logodata['stacklabel'].unique()
    if ncol is None:
        ncol = len(facets)

    # generate list of matrices to facet
    matrices = []
    for f in facets:
        facetdata = (logodata.query('facetlabel == @f')
                     .drop('facetlabel', axis=1)
                     .set_index('stacklabel')
                     .reindex(stacks)
                     .fillna(0)
                     )
        m = r.matrix(
                facetdata.values.ravel(),
                ncol=len(stacks),
                dimnames=[letters, stacks]
                )
        matrices.append(m)
    matrices = ListVector(TaggedList(matrices,
            tags=facets.astype('str')))

    # define the plotting function in a string
    rfuncstr = """
        facetedGGSeqLogo <- function(
            matrices, plotfile,
            ncol, width, height, 
            xname, xlabels, xlabelsrotate, xline, yname,
            letters, letter_colors
            )
        {
            p <- ggseqlogo(matrices, method='custom', ncol=ncol,
                    col_scheme=make_col_scheme(chars=letters,
                        cols=letter_colors)
                    ) +
                scale_x_continuous(xname, breaks=1:length(xlabels),
                    labels=xlabels) 

            if (xlabelsrotate) {
                axis.text.x = element_text(angle=90, hjust=1)
            } else {
                axis.text.x = element_text()
            }

            if (xline) {
                axis.line.x <- element_line(color='black')
            } else {
                axis.line.x <- element_blank()
            }

            if (nchar(trimws(yname))) {
                p <- p + scale_y_continuous(yname) +
                axis.line.y = element_line(color='black')
                axis.text.y = element_text()
            } else {
                axis.text.y <- element_blank()
                axis.line.y <- element_blank()
            }

            p <- p + theme(axis.text.x=axis.text.x,
                           axis.line.x=axis.line.x,
                           axis.text.y=axis.text.y,
                           axis.line.y=axis.line.y,
                           axis.text=element_text(size=12),
                           strip.text=element_text(size=13),
                           axis.title=element_text(size=13),
                           panel.spacing=unit(1.75, 'lines')
                           )

            ggsave(plotfile, plot=p, width=width, height=height)
        }
        """

    rfuncs = SignatureTranslatedAnonymousPackage(rfuncstr, "rfuncstr")
    rfuncs.facetedGGSeqLogo(
            matrices=matrices,
            plotfile=plotfile,
            ncol=ncol,
            width=width,
            height=height,
            xname='',
            xlabels=stacks,
            xlabelsrotate=True,
            xline=True,
            yname='',
            letters=StrVector(letters),
            letter_colors=StrVector([letter_colors[x] for x in letters])
            )

    if not os.path.isfile(plotfile):
        raise RuntimeError("failed to create {0}".format(plotfile))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
