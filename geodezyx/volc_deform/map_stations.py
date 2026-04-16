#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Map of GNSS stations from xyz_dic (ECEF => lat/lon)
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
import contextily as ctx
from pyproj import Transformer
from geodezyx import conv

xyz_dic = {
    "BOMG": [3352122.640732463, 4915723.753034734, -2297076.907396815],
    "BORG": [3352179.468244755, 4915389.35801885, -2297785.106621451],
    # "C98G": [3349647.071016342, 4918290.539762681, -2292685.261361883],
    # "CASG": [3342033.023264375, 4922773.708655806, -2290074.971696334],
    "CFNG": [3349537.184260634, 4916753.630176714, -2296761.787795528],
    # "CRAG": [3347500.97944078, 4917914.352085327, -2294826.506682635],
    "DERG": [3351077.500251723, 4916251.558666458, -2297303.247570688],
    "DSRG": [3351348.575051963, 4915747.889372488, -2298064.221329384],
    # "ENCG": [3354597.116872568, 4913638.028335221, -2297376.880814932],
    # "ENOG": [3353268.038184023, 4914797.26627217, -2296549.933107008],
    # "FEUG": [3340804.842314329, 4916311.646266831, -2305336.814085809],
    # "FJAG": [3350995.968928102, 4916475.044280337, -2295856.459849963],
    # "FOAG": [3350528.01283789, 4915139.525874629, -2299179.008076069],
    # "FREG": [3353758.5735742, 4915576.557677006, -2292887.943071418],
    # "GBNG": [3343297.13748501, 4919838.611709842, -2294437.898300031],
    # "GBSG": [3344238.331729389, 4917088.44551912, -2299848.399194048],
    "GITG": [3354343.97136548, 4914925.741250496, -2294766.595995225],
    # "GPNG": [3347812.884805933, 4917375.372526649, -2296671.256797257],
    # "GPSG": [3346186.407024792, 4916676.421682362, -2299378.496052825],
    # "HDLG": [3343663.622590312, 4918310.037725286, -2297365.60844327],
    # "MAIG": [3383479.378135485, 4901531.448891763, -2280412.059167621],
    # "PBRG": [3357258.662478508, 4913169.130239947, -2294728.333650563],
    # "PRAG": [3350809.181810517, 4913568.597968383, -2302266.991206395],
    # "PVDG": [3348379.438174965, 4916189.740796518, -2299156.356559711],
    # "RVAG": [3352538.751430386, 4914422.556831216, -2298369.016684549],
    "SNEG": [3351359.653965392, 4916209.867089082, -2296999.327082343],
    # "TRCG": [3342109.774250217, 4917359.715981821, -2301525.627671301],
}


def plot_stations_map(
    xyz_dic,
    zoom=13,
    margin_m=3000,
    figsize=None,
    dpi_max=150,
    dpi_min=72,
    marker_size=80,
    label_fontsize=7,
    cmap="terrain",
    tile_source=None,
    title="",
    save_path="/tmp/ovpf_stations_map.png",
):
    """
    Plot GNSS stations from an ECEF xyz dictionary on a topographic map.

    Parameters
    ----------
    xyz_dic : dict
        {site_name: [X, Y, Z]}  ECEF coordinates in metres.
    zoom : int
        OpenTopoMap tile zoom level (detail).
        Recommended range: 12 (broad) – 14 (fine).
    margin_m : int
        Map margin around the station bounding box, in metres (EPSG:3857).
        Increase to show more context / reduce apparent detail.
    figsize : tuple, int/float, or None
        Matplotlib figure size in inches.
        - ``None``        → width=10, height auto from geographic aspect ratio.
        - ``int / float`` → sets the **width**; height is derived automatically.
        - ``(w, h)``      → explicit width and height, no auto-adjustment.
    dpi_max : int
        DPI ceiling. The actual DPI is automatically lowered if needed so that
        native tile resolution > figure pixel resolution (prevents aliasing).
    dpi_min : int
        DPI floor (default 72). Prevents the auto-DPI from dropping so low that
        markers become invisible (happens with small margin_m / low zoom).
    marker_size : int or float
        Scatter marker area in points² (matplotlib ``s`` parameter). Default 80.
    label_fontsize : int or float
        Font size of station name labels in points. Default 7.
    cmap : str
        Matplotlib colormap for station height colouring.
    tile_source : contextily provider or None
        Tile source. Defaults to ctx.providers.OpenTopoMap.
    title : str
        Figure title.
    save_path : str or None
        Output file path. Set to None to skip saving.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    if tile_source is None:
        tile_source = ctx.providers.OpenTopoMap

    # --- Convert ECEF → lat/lon/h ---
    t_fwd = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    t_inv = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)

    sites, lats, lons, heigs = [], [], [], []
    for name, (x, y, z) in xyz_dic.items():
        lat, lon, h = conv.xyz2geo(x, y, z, outdeg=True)
        sites.append(name)
        lats.append(lat)
        lons.append(lon)
        heigs.append(h)

    lats  = np.array(lats)
    lons  = np.array(lons)
    heigs = np.array(heigs)

    # --- Reproject to Web Mercator (native tile CRS → zero reprojection aliasing) ---
    xs, ys = t_fwd.transform(lons, lats)

    x_min, x_max = xs.min() - margin_m, xs.max() + margin_m
    y_min, y_max = ys.min() - margin_m, ys.max() + margin_m

    # --- Auto figsize: preserve geographic aspect ratio (compensates for margins) ---
    # Accept figsize=None, figsize=<width>, or figsize=(w, h)
    if isinstance(figsize, (int, float)):
        fig_w = float(figsize)
        figsize = None          # trigger auto-height below
    else:
        fig_w = 10.0            # default width when figsize is None

    if figsize is None:
        aspect  = (y_max - y_min) / (x_max - x_min)   # data aspect (Web Mercator)
        # ~80 % of fig width goes to the axes (rest = colorbar + margins)
        axes_w  = fig_w * 0.80
        axes_h  = axes_w * aspect
        # ~88 % of fig height goes to the axes (rest = title + tick labels)
        fig_h   = round(axes_h / 0.88, 2)
        figsize = (fig_w, fig_h)
        print(f"auto figsize={figsize}  (data aspect={aspect:.3f})")

    # --- DPI auto-adjusted so native tile px ≥ figure px (always downsampling) ---
    n_tiles_x  = math.ceil((x_max - x_min) / (40075016.68 / 2**zoom))
    native_px  = n_tiles_x * 256
    dpi        = max(dpi_min, min(dpi_max, int(0.95 * native_px / figsize[0])))
    print(f"zoom={zoom}  native_px={native_px}  dpi_used={dpi}  "
          f"output={int(figsize[0]*dpi)}×{int(figsize[1]*dpi)} px")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.set_aspect("equal")   # 1 m = 1 m in x and y → no map distortion

    sc = ax.scatter(xs, ys, c=heigs, cmap=cmap, s=marker_size,
                    zorder=5, edgecolors="k", linewidths=0.5)

    label_offset = (x_max - x_min) * 0.01   # 1 % of the view extent
    for name, x, y in zip(sites, xs, ys):
        ax.text(x + label_offset, y + label_offset, name,
                fontsize=label_fontsize, zorder=6, ha="left", va="bottom",
                bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.6, lw=0))

    ctx.add_basemap(ax, source=tile_source, zoom=zoom,
                    interpolation="sinc", zorder=0, reset_extent=False)

    def fmt_lon(x, pos):
        lon, _ = t_inv.transform(x, 0)
        return f"{lon:.3f}°"

    def fmt_lat(y, pos):
        _, lat = t_inv.transform(0, y)
        return f"{lat:.3f}°"

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(fmt_lon))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(fmt_lat))
    ax.grid(True, linestyle="--", linewidth=0.4, color="grey", alpha=0.6, zorder=1)

    cbar = plt.colorbar(sc, ax=ax, orientation="vertical", pad=0.02, shrink=0.85)
    cbar.set_label("Height (m)")

    ax.set_title(title, fontsize=12)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig, ax


def plot_baselines_map(
    df_bl,
    xyz_dic,
    col="strain",
    agg="last",
    n_last=1,
    zoom=13,
    margin_m=3000,
    figsize=None,
    dpi_max=150,
    dpi_min=72,
    marker_size=80,
    label_fontsize=7,
    cmap_stations="terrain",
    tile_source=None,
    title="",
    save_path=None,
    cmap_lines="RdBu_r",
    line_width=4,
    line_alpha=0.85,
    show_values=True,
    value_fmt="{:.2e}",
    figax_tup=None,
):
    """
    Plot GNSS baselines as coloured line segments on a topographic map.

    Internally calls :func:`plot_stations_map` to draw the background map and
    station markers, then overlays one coloured segment per baseline pair.

    Parameters
    ----------
    df_bl : pandas.DataFrame
        Baselines DataFrame as returned by :func:`calc_baselines_direct` or
        :func:`calc_baselines_virtual`.  Must contain at least the columns
        ``site1``, ``site2``, and the column named by *col*.
    xyz_dic : dict
        ``{site_name: [X, Y, Z]}``  ECEF coordinates in metres.
    col : str
        Column of *df_bl* to map onto line colour.  Default ``"strain"``.
    agg : str or callable
        How to aggregate the per-epoch values of *col* into one scalar per
        baseline.  Supported strings: ``"last"``, ``"mean"``, ``"median"``.
        Any callable ``f(pd.Series) -> float`` is also accepted.
        Default is ``"last"``.
    n_last : int
        When *agg* is ``"last"``, average the last *n_last* non-NaN samples
        instead of the single latest value.  Ignored for other aggregators.
        Default is 1.
    zoom : int
        OpenTopoMap tile zoom level forwarded to :func:`plot_stations_map`.
    margin_m : int
        Map margin in metres around the station bounding box.
    figsize : tuple, int/float, or None
        Figure size in inches (forwarded to :func:`plot_stations_map`).
    dpi_max, dpi_min : int
        DPI ceiling and floor (forwarded to :func:`plot_stations_map`).
    marker_size : int or float
        Scatter marker area in points² for station markers.
    label_fontsize : int or float
        Font size for station name labels and optional midpoint value labels.
    cmap_stations : str
        Colormap for station height colouring.
    tile_source : contextily provider or None
        Tile source for the basemap.
    title : str
        Figure title.
    save_path : str or None
        Output file path.  ``None`` skips saving.
    cmap_lines : str
        Diverging colormap used to colour the baseline segments.
        Default ``"RdBu_r"``.
    line_width : int or float
        Line width of the baseline segments.  Default 4.
    line_alpha : float
        Opacity of the baseline segments.  Default 0.85.
    show_values : bool
        If ``True``, annotate the midpoint of each segment with the aggregated
        scalar value.  Default ``True``.
    value_fmt : str
        Python format string used when *show_values* is ``True``.
        Default ``"{:.2e}"``.
    figax_tup : tuple or None
        Existing ``(fig, ax)`` produced by a previous call to
        :func:`plot_stations_map`.  When provided the basemap step is skipped
        and lines are added directly.  Default ``None``.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    # ── 1. Aggregate df_bl per (site1, site2) → scalar ──────────────────────
    scalars = {}
    for (s1, s2), grp in df_bl.groupby(["site1", "site2"]):
        ser = grp[col].dropna()
        if len(ser) == 0:
            scalars[(s1, s2)] = np.nan
            continue
        if agg == "last":
            scalars[(s1, s2)] = float(ser.iloc[-n_last:].mean())
        elif agg == "mean":
            scalars[(s1, s2)] = float(ser.mean())
        elif agg == "median":
            scalars[(s1, s2)] = float(ser.median())
        elif callable(agg):
            scalars[(s1, s2)] = float(agg(ser))
        else:
            scalars[(s1, s2)] = float(ser.agg(agg))

    # ── 2. Build / reuse basemap figure ─────────────────────────────────────
    if figax_tup is None:
        fig, ax = plot_stations_map(
            xyz_dic,
            zoom=zoom,
            margin_m=margin_m,
            figsize=figsize,
            dpi_max=dpi_max,
            dpi_min=dpi_min,
            marker_size=marker_size,
            label_fontsize=label_fontsize,
            cmap=cmap_stations,
            tile_source=tile_source,
            title=title,
            save_path=None,   # save only after adding lines
        )
    else:
        fig, ax = figax_tup

    # ── 3. Reproject stations to Web Mercator ────────────────────────────────
    t_fwd = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    xy_merc = {}
    for name, (x, y, z) in xyz_dic.items():
        lat, lon, h = conv.xyz2geo(x, y, z, outdeg=True)
        xm, ym = t_fwd.transform(lon, lat)
        xy_merc[name] = (xm, ym)

    # ── 4. Colour normalization (symmetric around 0) ─────────────────────────
    valid_vals = np.array([v for v in scalars.values() if not np.isnan(v)])
    if len(valid_vals) == 0:
        print("plot_baselines_map: no valid scalar values — returning empty map.")
        return fig, ax

    vmax = np.abs(valid_vals).max()
    if vmax == 0:
        norm = mcolors.Normalize(vmin=-1, vmax=1)
    else:
        norm = mcolors.TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    cmap_obj = mcm.get_cmap(cmap_lines)

    # ── 5. Draw coloured line segments ───────────────────────────────────────
    view_w = ax.get_xlim()[1] - ax.get_xlim()[0]
    label_offset = view_w * 0.01

    for (s1, s2), val in scalars.items():
        if s1 not in xy_merc or s2 not in xy_merc:
            continue
        if np.isnan(val):
            continue
        x1, y1 = xy_merc[s1]
        x2, y2 = xy_merc[s2]
        color = cmap_obj(norm(val))
        ax.plot(
            [x1, x2], [y1, y2],
            color=color, lw=line_width, alpha=line_alpha,
            zorder=4, solid_capstyle="round",
        )
        if show_values:
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2
            ax.text(
                mx, my, value_fmt.format(val),
                fontsize=label_fontsize, zorder=7,
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.75, lw=0),
            )

    # ── 6. Colorbar for baseline metric ──────────────────────────────────────
    sm = mcm.ScalarMappable(cmap=cmap_obj, norm=norm)
    sm.set_array([])
    cbar2 = plt.colorbar(sm, ax=ax, orientation="vertical",
                         pad=0.08, shrink=0.7, anchor=(0.0, 0.35))
    cbar2.set_label(col)

    if save_path:
        plt.savefig(save_path, dpi=fig.dpi, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig, ax


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    plot_stations_map(
        xyz_dic,
        zoom=14,
        dpi_max=300,
        margin_m=400,
        save_path="/tmp/ovpf_stations_map.png",
        marker_size=150,
        label_fontsize=15,
        figsize=15,
        title="GNSS stations – Piton de la Fournaise (OVPF)"
    )
    plt.show()
