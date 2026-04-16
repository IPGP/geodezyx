import pandas as pd
import matplotlib.pyplot as plt

from geodezyx import utils, conv, files_rw, time_series, reffram, stats, utils_xtra
import numpy as np


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


def calc_baselines_virtual(df_inp, rov12_pairs, pivots, threshold_mad=3.5,
                           mean_win=86400, strain_win=7 * 86400):
    """
    Compute virtual baselines between rover pairs via a common pivot station.

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
        Concatenated DataFrame of baseline metrics (as returned by
        :func:`_d_calc`) with additional columns:

        * ``epoch``  : Observation epoch.
        * ``site1``  : First rover site name.
        * ``site2``  : Second rover site name.
        * ``pivot``  : Pivot station used to derive the virtual baseline.
    """

    pivots = utils.listify(pivots)
    rov12_sets_use = [set(e) for e in rov12_pairs]
    rov12_sets_done = []
    ii = -1
    col = ["x", "y", "z"]
    df_bl_stk = []

    df_pivs = df_inp[df_inp["base"].isin(pivots)]

    for piv, df_piv in df_pivs.groupby("base"):
        for rov1, df_rov1_ini in df_piv.groupby("rover"):
            for rov2, df_rov2_ini in df_piv.groupby("rover"):
                rov12_set = {rov1, rov2}
                if rov1 == rov2:
                    continue
                if not rov12_set in rov12_sets_use:
                    continue
                if rov12_set in rov12_sets_done:
                    continue

                df_rov1_wrk = df_rov1_ini.drop_duplicates(keep="first")
                df_rov2_wrk = df_rov2_ini.drop_duplicates(keep="first")

                thd = threshold_mad
                df_rov1_wrk, _ = stats.outlier_mad_df(df_rov1_wrk, col, thd)
                df_rov2_wrk, _ = stats.outlier_mad_df(df_rov2_wrk, col, thd)

                df_rov1_cmn, df_rov2_cmn = reffram.orb_df_common_epoch_finder(
                    df_rov1_wrk, df_rov2_wrk, order=["epoch"]
                )

                coords_delta = df_rov1_cmn[col] - df_rov2_cmn[col]
                df_bl = d_calc(coords_delta, mean_win=mean_win, strain_win=strain_win)

                df_bl["epoch"] = df_rov1_cmn.index
                df_bl["site1"] = rov1
                df_bl["site2"] = rov2
                df_bl["pivot"] = piv

                df_bl, _ = stats.outlier_mad_df(df_bl, ["d"], thd)

                df_bl_stk.append(df_bl)
                rov12_sets_done.append(rov12_set)

    df_bls = pd.concat(df_bl_stk)
    return df_bls


def calc_baselines_direct(
    df_inp, rovbas_pairs, bases_excluded=[], threshold_mad=3.5, xyz_dic_inp=None, mean_win=86400, strain_win=7 * 86400
):
    """
    Compute direct baselines between rovers and their reference base stations.

    Each rover position is differenced from the known (or first-epoch) position
    of its base station to obtain a baseline displacement time series.

    Parameters
    ----------
    df_inp : pandas.DataFrame
        Input DataFrame containing at least the columns ``rover``, ``base``,
        ``epoch``, ``x``, ``y``, ``z``.
    rovbas_pairs : list of tuple
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
        Concatenated DataFrame of baseline metrics (as returned by
        :func:`_d_calc`) with additional columns:

        * ``epoch``  : Observation epoch.
        * ``site1``  : Rover site name.
        * ``site2``  : Base station name.
        * ``pivot``  : Always ``None`` for direct baselines.
    """

    df_bls = []
    ii = -1
    col = ["x", "y", "z"]
    df_bl_stk = []
    thd = threshold_mad

    for (rov, bas), df_roba in df_inp.groupby(["rover", "base"]):

        if (rov, bas) not in rovbas_pairs:
            continue

        if bas in bases_excluded:
            continue

        rovbas_tup = (rov, bas)
        rovbas_set = set(rovbas_tup)

        df_wrk = df_roba
        ii += 1

        df_wrk, _ = stats.outlier_mad_df(df_wrk, col, thd)

        if not xyz_dic_inp or not bas in xyz_dic_inp:
            warnmsg = f"no reference coords for {bas}, will use the first one as ref."
            log.warning(warnmsg)
            coords_delta = df_wrk[col] - df_wrk[col].iloc[0]
        else:
            coords_delta = df_wrk[col] - xyz_dic_inp[bas]
        df_bl = d_calc(coords_delta, mean_win=mean_win, strain_win=strain_win)
        df_bl["epoch"] = df_wrk["epoch"].values
        df_bl["site1"] = rov
        df_bl["site2"] = bas
        df_bl["pivot"] = None

        df_bl, _ = stats.outlier_mad_df(df_bl, ["d"], thd)

        df_bl_stk.append(df_bl)

    df_bls = pd.concat(df_bl_stk)
    return df_bls


def baselines_plot(
    df_bl_inp,
    col="d_mean0",
    figax_tup=None,
    marker="",
    linestyle="-",
    suptitle="Direct baselines",
    ylabel = "Distance difference (cm)",
    plt_shift=0.02,
    plt_factor=100,
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
            df_plot["epoch"],
            (df_plot[col] + ii * plt_shift) * plt_factor,
            label=f"{sit1}-{sit2}",
            color="C" + str(ii),
            marker=marker,
            linestyle=linestyle,
        )
        ii += 1
        ax.set_ylabel(ylabel)
        ax.legend()
        fig.suptitle(suptitle)

    return fig, ax
