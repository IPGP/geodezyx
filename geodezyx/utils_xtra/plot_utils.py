#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.utils contains functions for operations
related to Python's plot operations.

it can be imported directly with:
from geodezyx import utils

The geodezyx toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/IPGP/geodezyx
"""

#### Import the logger
import logging
import os

########## BEGIN IMPORT ##########
#### External modules
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#### geodeZYX modules
from geodezyx import utils

log = logging.getLogger("geodezyx")

##########  END IMPORT  ##########


def color_list(l, colormap="jet"):
    """
    Generate a list of colors from a colormap for each unique value in L.

    Parameters
    ----------
    l : array-like
        Input list or array of values. The number of unique values determines the number of colors.
    colormap : str, optional
        Name of the matplotlib colormap to use. The default is 'jet'.

    Returns
    -------
    colist : list of tuple
        List of RGBA color tuples from the specified colormap.

    Notes
    -----
    The colors are evenly distributed across the colormap based on the number of unique values in L.

    See Also
    --------
    matplotlib.pyplot.get_cmap : Get a colormap by name.
    colors_from_colormap_getter : Alternative function to get colors from a colormap.
    """
    cm = plt.get_cmap(colormap)
    ncol = len(np.unique(l))
    colist = [cm(1.0 * i / ncol) for i in range(ncol)]
    return colist


def symbols_list(l=None):

    lsym = [
        "o",
        "v",
        "^",
        "<",
        ">",
        ".",
        ",",
        "1",
        "2",
        "3",
        "4",
        "8",
        "s",
        "p",
        "p",
        "*",
        "h",
        "H",
        "+",
        "x",
        "X",
        "D",
        "d",
        "|",
        "_",
    ]

    if not l:
        return lsym
    else:
        return lsym[: len(l)]


def colors_from_colormap_getter(ncolors, colormap="viridis"):
    """
    Get a list of colors from a matplotlib colormap.

    Parameters
    ----------
    ncolors : int
        Number of colors to generate.
    colormap : str, optional
        Name of the matplotlib colormap to use. The default is 'viridis'.

    Returns
    -------
    list of tuple
        List of RGBA color tuples from the specified colormap.

    See Also
    --------
    color_list : Generate colors based on unique values in a list.
    """
    import matplotlib.pyplot as plt

    cm = plt.get_cmap(colormap)
    return [cm(1.0 * i / ncolors) for i in range(ncolors)]


def ylim_easy(lin, delta=0.1, min_null_if_neg=False):
    """
    Calculate convenient axis limits for a data array.

    Parameters
    ----------
    lin : array-like
        Input data.
    delta : float, optional
        Fraction of the data range to add as padding. The default is 0.1.
    min_null_if_neg : bool, optional
        If True, set the lower limit to 0 if it would be negative. The default is False.

    Returns
    -------
    tuple of float
        (lower_limit, upper_limit) for the y-axis.

    Notes
    -----
    Useful for automatic axis limit calculation in plots.
    """
    minn = np.min(lin)
    maxx = np.max(lin)
    rangee = np.abs(maxx - minn)
    if min_null_if_neg and (minn - delta * rangee) < 0.0:
        return (0.0, maxx + delta * rangee)
    else:
        return (minn - delta * rangee, maxx + delta * rangee)


def get_figure(figin=0):
    """
    Get or create a matplotlib Figure object.

    Parameters
    ----------
    figin : int or matplotlib.figure.Figure, optional
        Figure specification. If 0, creates a new figure.
        If an integer, returns figure with that number.
        If a Figure object, returns that figure. The default is 0.

    Returns
    -------
    figout : matplotlib.figure.Figure
        The requested or created figure object with at least one axes.

    Notes
    -----
    Ensures the returned figure has at least one axes (subplot).
    """
    if isinstance(figin, matplotlib.figure.Figure):
        figout = figin
    elif figin == 0:
        figout = plt.figure()
    else:
        figout = plt.figure(figin)
    # be sure the fig have a axe (necessary ?)
    if figout.get_axes() == []:
        figout.add_subplot(111)
    plt.figure(figout.number)
    return figout


def figure_saver(
    figobjt_in,
    outdir,
    outname,
    outtype=(".png", ".pdf", ".figpik"),
    formt=None,
    dpi=200,
    transparent=False,
):
    """
    This function provides a front end to export pretty-print plots

    Parameters
    ----------
    figobjt_in : matplotlib Figure object
        input matplotlib Figure object. use for instance plt.gcf() to get it.
    outdir : str
        the output directory.
    outname : str
        output prefix filename.
    outtype : tuple, optional
        the output formats. The default is ('.png','.pdf','.figpik').
    formt : 2-tuple or string , optional
        the format (size) of the plot.
        if string: a Ax format (A4, A3 etc...)
        if tuple: size of the plot  in inches.
        The default is None.
    dpi : int, optional
        DPI of the figure. The default is 200.
    transparent : bool, optional
        make the plot transparent. The default is False.

    Returns
    -------
    outpath_stk : string or list of string
        output paths of the plots.

    """

    if not utils.is_iterable(outtype):
        outtype = (outtype,)

    outpath_stk = []
    for outtype_iter in outtype:
        if "pik" in outtype_iter:
            outpath = utils.pickle_saver(figobjt_in, outdir, outname, outtype_iter)
        else:
            outpath = os.path.join(outdir, outname + outtype_iter)

            formtup = None
            if formt:
                if type(formt) is tuple:
                    formtup = formt
                elif type(formt) is str:
                    if formt.upper() == "A4":
                        formtup = (11.69, 8.27)
                    elif formt.upper() == "A3":
                        formtup = (16.53, 11.69)
                    else:
                        log.warning("assume Figure format as A4")
                        formtup = (11.69, 8.27)

                if formtup:
                    figobjt_in.set_size_inches(*formtup)

            figobjt_in.savefig(outpath, transparent=transparent, dpi=dpi)

        outpath_stk.append(outpath)

    if len(outpath_stk) == 1:
        outpath_stk = outpath_stk[0]
    return outpath_stk


def axis_data_coords_sys_transform(axis_obj_in, xin, yin, inverse=False):
    """inverse = False : Axis => Data
    = True  : Data => Axis
    """
    xlim = axis_obj_in.get_xlim()
    ylim = axis_obj_in.get_ylim()

    xdelta = xlim[1] - xlim[0]
    ydelta = ylim[1] - ylim[0]
    if not inverse:
        xout = xlim[0] + xin * xdelta
        yout = ylim[0] + yin * ydelta
    else:
        xdelta2 = xin - xlim[0]
        ydelta2 = yin - ylim[0]
        xout = xdelta2 / xdelta
        yout = ydelta2 / ydelta
    return xout, yout


def id2val(value_lis, id_lis, idin):
    """from a value list and a id pointer list
    return the good val from the good id
    replace dico bc. set is not supproted as key"""
    return value_lis[id_lis.index(idin)]


def set_size_for_pub(width=418.25368, fraction=1, subplot=None):
    """
    Set aesthetic figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width : float, optional
        Width of the figure in points (pt). The default is 418.25368 (approximately 146 mm).
    fraction : float, optional
        Fraction of the width that the figure should occupy. The default is 1.
    subplot : list of int, optional
        Subplot grid dimensions as [nrows, ncols]. The default is [1, 1].

    Returns
    -------
    fig_dim : tuple of float
        Dimensions of the figure as (width_inches, height_inches).

    Notes
    -----
    Uses the golden ratio (φ = (√5 - 1) / 2 ≈ 0.618) to set aesthetic figure height.
    Useful for creating publication-ready plots that fit nicely in LaTeX documents.

    Examples
    --------
    >>> width = 418.25368  # Standard LaTeX column width
    >>> fig_dim = set_size_for_pub(width, fraction=0.5, subplot=[2, 2])
    >>> fig = plt.figure(figsize=fig_dim)
    """
    if subplot is None:
        subplot = [1, 1]
    # Width of figure
    fig_width_pt = width * fraction

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5**0.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplot[0] / subplot[1])

    fig_dim = (fig_width_in, fig_height_in)

    return fig_dim


def gaussian_for_plot(d, density=False, nbins=500, nsigma=3.5):
    """
    Generate a Gaussian curve for histogram overlay plots.

    Parameters
    ----------
    d : array-like
        Data vector to fit a Gaussian distribution to.
    density : bool, optional
        If True, returns the PDF (normalized). If False, scales the PDF to match histogram area.
        The default is False.
    nbins : int, optional
        Number of bins (points) to generate for the curve. The default is 500.
    nsigma : float, optional
        Number of standard deviations to span for the x-axis (μ ± nsigma*σ).
        The default is 3.5.

    Returns
    -------
    xpdf : numpy.ndarray
        X coordinates of the Gaussian curve.
    ypdf_out : numpy.ndarray
        Y coordinates of the Gaussian curve (PDF or histogram-scaled).

    Notes
    -----
    Useful for overlaying a fitted Gaussian curve on a histogram.
    When density=False, the curve is scaled to match the histogram's area.
    When density=True, the curve is the probability density function.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> data = np.random.randn(1000)
    >>> x_curve, y_curve = gaussian_for_plot(data)
    >>> plt.hist(data, bins=50, density=True)
    >>> plt.plot(x_curve, y_curve, 'r-', label='Gaussian fit')

    See Also
    --------
    scipy.stats.norm.pdf : Probability density function for normal distribution.
    """
    import scipy
    from scipy import integrate

    mu = np.mean(d)
    sigma = np.std(d)
    xpdf = np.linspace(mu - nsigma * sigma, mu + nsigma * sigma, nbins)
    ypdf = scipy.stats.norm.pdf(xpdf, mu, sigma)
    ypdf_out = ypdf
    if not density:
        ybin, xbin = np.histogram(d, bins=nbins)
        area_bin = integrate.trapezoid(ybin, dx=np.diff(xbin)[0])

        ypdf_out = ypdf * area_bin
    return xpdf, ypdf_out
