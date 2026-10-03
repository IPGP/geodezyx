import pandas as pd
import matplotlib.pyplot as plt

from geodezyx import utils, conv, files_rw, time_series, reffram, stats, utils_xtra
import numpy as np

#### Import the logger
import logging

log = logging.getLogger("geodezyx")


def d_calc(coords_delta, mean_win=86400, strain_win=7 * 86400):
    """
    Compute baseline distance metrics and strain rate from coordinate differences.

    Parameters
    ----------
    coords_delta : array-like of shape (N, 3)
        Array of coordinate differences (dx, dy, dz) for each epoch.
    mean_win : int, optional
        Rolling window size (in samples) used to compute the running mean of
        the baseline distance. Default is 86400 (e.g. 1 day at 1 Hz).
    strain_win : int or list of int, optional
        Window size(s) (in samples) used to compute the strain rate as the
        relative change of the running mean distance over the window.
        Default is 7 * 86400 (e.g. 7 days at 1 Hz).

    Returns
    -------
    df_bl : pandas.DataFrame
        DataFrame with the following columns:

        * ``d``        : Euclidean baseline distance.
        * ``d0``       : Distance centred on its median.
        * ``d_mean``   : Rolling mean of the distance.
        * ``d_mean0``  : Rolling mean centred on the overall median.
        * ``d_diff``   : First difference of the distance.
        * ``strain``   : Strain rate computed over the first window in
                         *strain_win* (relative change: (last - first) / first).
        * ``strain<w>`` : (only when multiple windows are given) Strain rate
                          for each additional window size *w*.
    """

    strain_win = utils.listify(strain_win)

    d = np.linalg.norm(coords_delta, axis=1)
    d_med = np.median(d)
    d_ser = pd.Series(d)

    d_mean = d_ser.rolling(window=mean_win).mean().values
    d_diff = np.diff(d, prepend=d[0])  # or append=d[-1]

    df_bl = pd.DataFrame()

    df_bl["d"] = d
    df_bl["d0"] = d - d_med
    df_bl["d_mean"] = d_mean
    df_bl["d_mean0"] = d_mean - d_med
    df_bl["d_diff"] = d_diff

    # Compute linear regression
    a, b = stats.linear_regression(df_bl["d_mean"].values, df_bl["epoch"].values)

    # strain: relative change of d_mean between first and last sample of the strain_win window
    lbd_strain = lambda x: (x[-1] - x[0]) / x[0]
    d_mean_ser = pd.Series(d_mean)

    for isw, sw in enumerate(strain_win):
        strain = d_mean_ser.rolling(window=sw).apply(lbd_strain, raw=True).values
        if isw == 0:
            df_bl["strain"] = strain
        if len(strain_win) > 1:
            df_bl["strain" + str(sw)] = strain

    return df_bl


def _process_bl_pair(df_bl, thd):
    """Apply common baseline filtering."""


    df_bl_cor = df_bl["d"] - df_bl["d_mean"]
    df_bl_cor = pd.DataFrame(df_bl_cor)
    df_bl_cor.columns = ["d"]
    _, mad_bool = stats.outlier_mad_df(df_bl_cor, ["d"], thd)
    return df_bl[mad_bool]


def _resample_bl(df_bl_inp, resample_inp):
    df_bl_inp.set_index("epoch", inplace=True)
    df_bl_inp = df_bl_inp.resample(resample_inp).median(numeric_only=True)
    df_bl_inp.reset_index(inplace=True)
    return df_bl_inp


def _virtual_baseline(
    df_inp, rov12_pairs, pivots, threshold_mad, mean_win, strain_win, resample
):
    """Compute virtual baselines."""
    col = ["x", "y", "z"]
    df_bl_stk = []

    pivots = utils.listify(pivots)
    rov12_sets_use = [set(e) for e in rov12_pairs]
    rov12_sets_done = []

    df_pivs = df_inp[df_inp["base"].isin(pivots)]

    for piv, df_piv in df_pivs.groupby("base"):
        for rov1, df_rov1_ini in df_piv.groupby("rover"):
            for rov2, df_rov2_ini in df_piv.groupby("rover"):
                rov12_set = {rov1, rov2}
                if rov1 == rov2 or rov12_set in rov12_sets_done:
                    continue
                if not rov12_set in rov12_sets_use:
                    continue

                df_rov1_wrk = _resample_bl(
                    df_rov1_ini.drop_duplicates(keep="first"), resample
                )
                df_rov2_wrk = _resample_bl(
                    df_rov2_ini.drop_duplicates(keep="first"), resample
                )

                df_rov1_wrk, _ = stats.outlier_mad_df(df_rov1_wrk, col, threshold_mad)
                df_rov2_wrk, _ = stats.outlier_mad_df(df_rov2_wrk, col, threshold_mad)

                df_rov1_cmn, df_rov2_cmn = reffram.orb_df_common_epoch_finder(
                    df_rov1_wrk, df_rov2_wrk, order=["epoch"]
                )

                coords_delta = df_rov1_cmn[col] - df_rov2_cmn[col]
                df_bl = d_calc(coords_delta, mean_win=mean_win, strain_win=strain_win)
                df_bl["epoch"] = df_rov1_cmn.index
                df_bl["site1"] = rov1
                df_bl["site2"] = rov2
                df_bl["pivot"] = piv

                df_bl = _process_bl_pair(df_bl, threshold_mad)
                df_bl_stk.append(df_bl)
                rov12_sets_done.append(rov12_set)

    return df_bl_stk


def _direct_baseline(
    df_inp,
    rov12_pairs,
    bases_excluded,
    threshold_mad,
    xyz_dic_inp,
    mean_win,
    strain_win,
    resample,
):
    """Compute direct baselines."""
    col = ["x", "y", "z"]
    df_bl_stk = []

    for (rov, bas), df_roba in df_inp.groupby(["rover", "base"]):
        if (rov, bas) not in rov12_pairs or bas in bases_excluded:
            continue

        df_wrk = _resample_bl(df_roba, resample)
        df_wrk, _ = stats.outlier_mad_df(df_wrk, col, threshold_mad)

        if not xyz_dic_inp or bas not in xyz_dic_inp:
            log.warning(
                f"no reference coords for {bas}, will use the first one as ref."
            )
            coords_delta = df_wrk[col] - df_wrk[col].iloc[0]
        else:
            coords_delta = df_wrk[col] - xyz_dic_inp[bas]

        # Compute baseline metrics
        df_bl = d_calc(coords_delta, mean_win=mean_win, strain_win=strain_win)
        # Add metadata columns
        df_bl["epoch"] = df_wrk["epoch"].values
        df_bl["site1"] = rov
        df_bl["site2"] = bas
        df_bl["pivot"] = None

        df_bl = _process_bl_pair(df_bl, threshold_mad)
        df_bl_stk.append(df_bl)

    return df_bl_stk


def calc_baselines(
    df_inp,
    mode="direct",
    rov12_pairs=None,
    pivots=None,
    bases_excluded=None,
    threshold_mad=3.5,
    xyz_dic_inp=None,
    mean_win=86400,
    strain_win=7 * 86400,
    resample="1min",
):
    """
    Compute baselines between rovers and reference stations.

    Can compute either direct baselines (rover vs base station) or virtual
    baselines (rover pairs via a common pivot station).

    Parameters
    ----------
    df_inp : pandas.DataFrame
        Input DataFrame containing at least the columns ``base``, ``rover``,
        ``epoch``, ``x``, ``y``, ``z``.
    mode : str, optional
        Baseline computation mode. Either ``"direct"`` or ``"virtual"``.
        Default is ``"direct"``.
    rov12_pairs : list of tuple or list of set, optional
        Pairs to process:
        - For ``mode="virtual"``: Rover pairs for which virtual baselines should
          be computed, e.g. ``[("ROV1", "ROV2"), ...]``.
        - For ``mode="direct"``: List of ``(rover, base)`` pairs, e.g.
          ``[("ROV1", "BASE1"), ...]``.
        Required for both modes.
    pivots : str or list of str, optional
        (For ``mode="virtual"``) Name(s) of the pivot (base) station(s) to use.
        Required if ``mode="virtual"``.
    bases_excluded : list of str, optional
        (For ``mode="direct"``) Base station names to skip.
        Default is an empty list.
    threshold_mad : float, optional
        Median Absolute Deviation multiplier used as the outlier rejection
        threshold. Default is 3.5.
    xyz_dic_inp : dict or None, optional
        (For ``mode="direct"``) Dictionary mapping base station name to its
        reference coordinates ``{base: [x, y, z]}``. If ``None`` or the base
        is not found in the dictionary the first epoch of the rover time series
        is used as reference. Default is ``None``.
    mean_win : int, optional
        Rolling window size (in samples) used to compute the running mean of
        the baseline distance. Default is 86400.
    strain_win : int or list of int, optional
        Window size(s) (in samples) used to compute the strain rate.
        Default is 7 * 86400.
    resample : str, optional
        Resampling frequency. Default is ``"1min"``.

    Returns
    -------
    df_bls : pandas.DataFrame
        Concatenated DataFrame of baseline metrics (as returned by
        :func:`d_calc`) with additional columns:

        * ``epoch``  : Observation epoch.
        * ``site1``  : First site name (rover for direct, first rover for virtual).
        * ``site2``  : Second site name (base for direct, second rover for virtual).
        * ``pivot``  : Pivot station used (``None`` for direct baselines).
    """

    if mode not in ("direct", "virtual"):
        raise ValueError(f"mode must be 'direct' or 'virtual', got '{mode}'")

    if rov12_pairs is None:
        raise ValueError("rov12_pairs is required")

    if bases_excluded is None:
        bases_excluded = []

    # Dispatch to mode-specific function
    if mode == "virtual":
        if pivots is None:
            raise ValueError("pivots is required when mode='virtual'")
        df_bl_stk = _virtual_baseline(
            df_inp, rov12_pairs, pivots, threshold_mad, mean_win, strain_win, resample
        )
    else:  # mode == "direct"
        df_bl_stk = _direct_baseline(
            df_inp,
            rov12_pairs,
            bases_excluded,
            threshold_mad,
            xyz_dic_inp,
            mean_win,
            strain_win,
            resample,
        )

    return pd.concat(df_bl_stk)


# Backward compatibility wrappers
def calc_baselines_virtual(
    df_inp,
    rov12_pairs,
    pivots,
    threshold_mad=3.5,
    mean_win=86400,
    strain_win=7 * 86400,
    resample="1min",
):
    """
    Compute virtual baselines between rover pairs via a common pivot station.

    .. deprecated::
        Use :func:`calc_baselines` with ``mode="virtual"`` instead.

    For each pivot station the positions of two rovers are differenced to
    produce a virtual baseline, removing common-mode errors introduced by the
    pivot.

    Parameters
    ----------
    df_inp : pandas.DataFrame
        Input DataFrame containing at least the columns ``base``, ``rover``,
        ``epoch``, ``x``, ``y``, ``z``.
    rov12_pairs : list of tuple or list of set
        Rover pairs for which virtual baselines should be computed, e.g.
        ``[("ROV1", "ROV2"), ...]``.
    pivots : str or list of str
        Name(s) of the pivot (base) station(s) to use.
    threshold_mad : float, optional
        Median Absolute Deviation multiplier used as the outlier rejection
        threshold. Default is 3.5.

    Returns
    -------
    df_bls : pandas.DataFrame
        Concatenated DataFrame of baseline metrics.
    """
    return calc_baselines(
        df_inp,
        mode="virtual",
        rov12_pairs=rov12_pairs,
        pivots=pivots,
        threshold_mad=threshold_mad,
        mean_win=mean_win,
        strain_win=strain_win,
        resample=resample,
    )


def calc_baselines_direct(
    df_inp,
    rov12_pairs,
    bases_excluded=None,
    threshold_mad=3.5,
    xyz_dic_inp=None,
    mean_win=86400,
    strain_win=7 * 86400,
    resample="1min",
):
    """
    Compute direct baselines between rovers and their reference base stations.

    .. deprecated::
        Use :func:`calc_baselines` with ``mode="direct"`` instead.

    Each rover position is differenced from the known (or first-epoch) position
    of its base station to obtain a baseline displacement time series.

    Parameters
    ----------
    df_inp : pandas.DataFrame
        Input DataFrame containing at least the columns ``rover``, ``base``,
        ``epoch``, ``x``, ``y``, ``z``.
    rov12_pairs : list of tuple
        List of ``(rover, base)`` pairs to process, e.g.
        ``[("ROV1", "BASE1"), ...]``.
    bases_excluded : list of str, optional
        Base station names to skip. Default is an empty list.
    threshold_mad : float, optional
        Median Absolute Deviation multiplier used as the outlier rejection
        threshold. Default is 3.5.
    xyz_dic_inp : dict or None, optional
        Dictionary mapping base station name to its reference coordinates
        ``{base: [x, y, z]}``. If ``None`` or the base is not found in the
        dictionary the first epoch of the rover time series is used as
        reference. Default is ``None``.

    Returns
    -------
    df_bls : pandas.DataFrame
        Concatenated DataFrame of baseline metrics.
    """
    if bases_excluded is None:
        bases_excluded = []
    return calc_baselines(
        df_inp,
        mode="direct",
        rov12_pairs=rov12_pairs,
        bases_excluded=bases_excluded,
        threshold_mad=threshold_mad,
        xyz_dic_inp=xyz_dic_inp,
        mean_win=mean_win,
        strain_win=strain_win,
        resample=resample,
    )


def baselines_plot(
    df_bl_inp,
    col="d_mean0",
    figax_tup=None,
    marker="",
    linestyle="-",
    suptitle="Direct baselines",
    ylabel="Distance difference (cm)",
    plt_shift=0.02,
    plt_factor=100,
    decim=100,
):
    """
    Plot baseline distance time series for all site pairs.

    Parameters
    ----------
    df_bl_inp : pandas.DataFrame
        DataFrame of baselines as returned by :func:`calc_baselines_direct` or
        :func:`calc_baselines_virtual`. Must contain the columns ``site1``,
        ``site2``, ``epoch`` and the column specified by *d_col*.
    col : str, optional
        Name of the column to plot on the y-axis. Default is ``"d_mean0"``
        (centred running mean distance).
    figax_tup : tuple of (matplotlib.figure.Figure, matplotlib.axes.Axes) or None, optional
        Existing ``(fig, ax)`` tuple to draw on. If ``None`` a new figure and
        axes are created. Default is ``None``.
    marker : str, optional
        Matplotlib marker style passed to ``ax.plot``. Default is ``""``
        (no marker).
    linestyle : str, optional
        Matplotlib line style passed to ``ax.plot``. Default is ``"-"``.
    suptitle : str, optional
        Title string for the figure. Default is ``"Direct baselines"``.
    ylabel : str, optional
        Label for the y-axis. Default is ``"Distance difference (cm)"``.
    plt_shift : float, optional
        Vertical shift applied to each site pair's time series to avoid overlap.
        Default is ``0.02`` (in the same units as the y-axis).
    plt_factor : float, optional
        Scaling factor applied to the y-axis values. Default is ``100``
        (to convert meters to centimeters).
    decim : int, optional
        Decimation factor for plotting. Only every *decim*-th point is plotted.
        Default is ``100``.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plot.
    ax : matplotlib.axes.Axes
        The axes object containing the plot.
    """

    fig, ax = figax_tup if figax_tup else plt.subplots()

    ii = 0
    for (sit1, sit2), df_plot in df_bl_inp.groupby(["site1", "site2"]):

        ax.plot(
            df_plot["epoch"].values[::decim],
            (df_plot[col].values[::decim] + ii * plt_shift) * plt_factor,
            label=f"{sit1}-{sit2}",
            color="C" + str(ii),
            marker=marker,
            linestyle=linestyle,
        )
        ii += 1

    last_epoc = df_bl_inp["epoch"].max()
    last_epoc_str = conv.dt2str_iso(last_epoc)
    now_str = conv.dt2str_iso(conv.now("utc"))
    ax.set_ylabel(ylabel)
    ax.legend()
    ax.set_title(f"generated: {now_str}, last epoch: {last_epoc_str}")
    fig.suptitle(suptitle)
    fig.tight_layout()

    return fig, ax
