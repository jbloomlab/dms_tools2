"""
===================
rplot
===================

Plotting that uses ``R``.

The `dms_tools2` software is written in Python, but there are
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
import warnings

import pandas

import phydmslib.weblogo
#: default colors for amino acid chars, by functional group
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
from rpy2.robjects.vectors import StrVector, ListVector, FloatVector
from rpy2.rlike.container import TaggedList
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

#: Show warnings when running ``R`` code? Plausible values: `default`,
#: `ignore`, and `always`. The `R` code gives many warnings, so
#: `ignore` is good when not developing new code.
SHOW_WARNINGS = 'ignore'

# install necessary packages
_packages = ['ggplot2', 'ggseqlogo']
utils = importr('utils')
with warnings.catch_warnings():
    warnings.simplefilter(SHOW_WARNINGS)
    utils.chooseCRANmirror(ind=1)
    utils.install_packages(StrVector(_packages))
for _package in _packages:
    importr(_package)

# read the R code that defines the R functions we run
_RCODEFILE = os.path.join(os.path.dirname(__file__), 'rplot_Rcode.R')
assert os.path.isfile(_RCODEFILE)
with open(_RCODEFILE) as f:
    _RCODE = f.read()
_RFUNCS = SignatureTranslatedAnonymousPackage(_RCODE, "_RFUNCS")


def versionInfo():
    """Returns string giving `rpy2` and ``R`` version info."""
    return "Using `rpy2` {0} with ``R`` {1}".format(
            rpy2.__version__, rinterface.R_VERSION_BUILD[1])


def siteSubsetGGSeqLogo(logodata, chars, plotfile, width, height,
        yname='', char_colors=AA_COLORS_FG, ylimits=None):
    """Creates one-row logo plot with subset of sites.

    Designed to show logo plot for a subset of sites. This
    is useful when you have data for many sites, but only
    want to look at a few of them. 

    Args:
        `logodata` (pandas DataFrame)
            Contains data to plot. Should have the columns
            `site`, `show`, and a column giving the height
            height of each char in `chars`. Only sites
            where `show` is `True` are shown. Sites are 
            shown in the order they occur in this dataframe,
            with spaces every time there is an interspersed
            site with `show` being `False`. 
        `chars` (list)
            Letters for which we plot heights.
        `plotfile` (str)
            Name of created plot.
        `width` (float)
            Width of plot in inches.
        `height` (float)
            Height of plot in inches.
        `yname` (str)
            If set to a non-empty string, is the y-axis label
            and yticks are drawn.
        `char_colors` (dict)
            Values give color for every character in `chars`.
        `ylimits` (`None` or 2-tuple)
            If not `None`, should give the ylimits for the plot
            as `(ymin, ymax)`

    Here is an example that creates a plot for a subset of
    sites for two characters:

    >>> logodata = pandas.read_csv(io.StringIO(
    ...     '''site show    A    C
    ...        A101 True  0.8  0.2
    ...        N102 True  0.7  0.3
    ...        K103 False 0.1  0.9
    ...        L104 True  0.8  0.2
    ...        S105 True  0.5  0.5
    ...        T106 False 0.2  0.8
    ...        G107 False 0.4  0.6
    ...        L108 True  0.7  0.3'''),
    ...     delim_whitespace=True, index_col=False)
    >>> plotfile = '_siteSubsetGGSeqLogo_test_plot.png'
    >>> siteSubsetGGSeqLogo(logodata,
    ...         chars=['A', 'C'],
    ...         plotfile=plotfile,
    ...         width=3.5, height=2
    ...         )
    >>> os.path.isfile(plotfile)
    True

    Here is the plot created by the code block above:

    .. image:: _static/_siteSubsetGGSeqLogo_test_plot.png
       :width: 55%
       :align: center

    """
    if os.path.isfile(plotfile):
        os.remove(plotfile)

    assert set(chars) <= set(char_colors.keys()), \
            "`char_colors` not defined for all chars"

    expectcol = ['site', 'show'] + chars
    assert set(logodata.columns) >= set(expectcol), \
            "`logodata` needs these column: {0}".format(expectcol)

    assert logodata['show'].any(), "no sites to show"

    # for each consecutive set of rows not to show, keep just one
    logodata = logodata[expectcol]
    logodata['keeprow'] = (
            ((logodata['show']) | 
                (logodata['show'] != logodata['show'].shift(1)))
            )
    logodata = logodata.query('keeprow').reset_index()

    # trim first and last row if they are not to be shown
    if not logodata.iloc[0]['show']:
        logodata = logodata.iloc[1 : ].reset_index()
    if not logodata.iloc[-1]['show']:
        logodata = logodata.iloc[ : -1]

    # set site label to empty and data to zero for rows not to show
    logodata.loc[~logodata['show'], 'site'] = ''
    logodata.loc[~logodata['show'], chars] = 0
    vertlines = logodata.query('~show').index.values + 1

    # generate matrix to plot
    sites = logodata['site']
    matrix = r.matrix(logodata.set_index('site')[chars].values.ravel(),
            ncol=len(sites),
            dimnames=[chars, sites]
            )

    if ylimits is None:
        ylimits = rinterface.NULL
    else:
        ylimits = FloatVector(ylimits)

    # make the plot
    with warnings.catch_warnings():
        warnings.simplefilter(SHOW_WARNINGS)
        _RFUNCS.siteSubsetGGSeqLogo(
                mat=matrix,
                plotfile=plotfile,
                width=width,
                height=height,
                xlabels=list(map(str, sites)),
                vertlines=vertlines,
                yname=yname,
                chars=StrVector(chars),
                char_colors=StrVector([char_colors[x] for x in chars]),
                ylimits=ylimits
                )

    if not os.path.isfile(plotfile):
        raise RuntimeError("failed to create {0}".format(plotfile))


def facetedGGSeqLogo(logodata, chars, plotfile, width, height,
        ncol=None, char_colors=AA_COLORS_FG, xlabelsrotate=True):
    """Creates faceted logo plot.

    Designed to show several measurements on the same site
    site-by-side, potentially for many sites. Each site
    must have the same set of measurements.

    Makes panel of logo plots faceted on `logodata['facetlabel']`,
    where character stacks are labeled by `logodata['stacklabel']`
    and show the characters at the indicated heights.

    Args:
        `logodata` (pandas DataFrame)
            Contains data to plot. Should have the columns
            `facetlabel`, `stacklabel`, and a column giving the
            height of each character in `chars`.
        `chars` (list)
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
        `char_colors` (dict)
            Values give color for every character in `chars`.
        `xlabelsrotate` (bool)
            Do we rotate the x-labels?

    Here is an example that creates two facets each with
    two stacks for the characters `A` and `C`:

    >>> logodata = pandas.read_csv(io.StringIO(
    ...     '''facetlabel  stacklabel   A   C
    ...            site-1       BF520 0.8 0.2
    ...            site-1       BG505 0.9 0.1
    ...            site-2       BF520 0.4 0.6
    ...            site-2       BG505 0.5 0.5'''),
    ...     delim_whitespace=True, index_col=False)
    >>> plotfile = '_facetedGGSeqLogo_test_plot.png'
    >>> facetedGGSeqLogo(logodata,
    ...         chars=['A', 'C'],
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

    assert set(chars) <= set(char_colors.keys()), \
            "`char_colors` not defined for all chars"

    # get and order data columns
    df_cols = ['facetlabel', 'stacklabel'] + chars
    assert set(logodata.columns) >= set(df_cols), "df lacks required columns"
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
                dimnames=[chars, stacks]
                )
        matrices.append(m)
    matrices = ListVector(TaggedList(matrices,
            tags=facets.astype('str')))

    # make the plot
    with warnings.catch_warnings():
        warnings.simplefilter(SHOW_WARNINGS)
        _RFUNCS.facetedGGSeqLogo(
                matrices=matrices,
                plotfile=plotfile,
                ncol=ncol,
                width=width,
                height=height,
                xname='',
                xlabels=stacks,
                xlabelsrotate=xlabelsrotate,
                xline=True,
                yname='',
                chars=StrVector(chars),
                char_colors=StrVector([char_colors[x] for x in chars])
                )

    if not os.path.isfile(plotfile):
        raise RuntimeError("failed to create {0}".format(plotfile))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
