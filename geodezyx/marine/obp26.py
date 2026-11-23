import numpy as np
import xarray as xr
import pandas as pd
import gsw


def calculate_bottom_pressure(
    ds,
    latitude,
    longitude,
    depth_sensor,
    temp_var,
    sal_var,
    ssh_var,
    temperature_type="potential",
    rho0=1025.0,
    g=9.80665,
):
    """
    Calculate non-tidal ocean bottom pressure anomalies.

    Parameters
    ----------
    ds : xarray.Dataset
        Ocean model dataset.

    latitude, longitude : float
        Station coordinates (degrees).

    depth_sensor : float
        Sensor depth below mean sea level (m).

    temp_var, sal_var, ssh_var : str
        Dataset variable names.

    temperature_type : str
        'potential' or 'insitu'.

    rho0 : float
        Reference seawater density (kg/m3).

    Returns
    -------
    xarray.Dataset
        Ocean bottom pressure anomalies.
    """

    # Select nearest model grid point
    lon = longitude

    if float(ds.longitude.max()) > 180 and lon < 0:
        lon = lon % 360

    point = ds.sel(
        latitude=latitude,
        longitude=lon,
        method="nearest",
    )

    # Extract coordinates
    z = point.depth.values.astype(float)

    # Ensure positive downward
    z = np.abs(z)

    order = np.argsort(z)
    z = z[order]

    # Extract temperature and salinity
    temp = point[temp_var].transpose("time", "depth").values[:, order]

    sal = point[sal_var].transpose("time", "depth").values[:, order]

    ssh = point[ssh_var].values.astype(float)

    # Check vertical coverage
    if depth_sensor > np.nanmax(z):
        raise ValueError("Sensor depth exceeds model vertical coverage")

    # Convert model depth to sea pressure (dbar)
    p = gsw.p_from_z(-z, latitude)

    # Practical Salinity -> Absolute Salinity
    SA = gsw.SA_from_SP(
        sal,
        p[np.newaxis, :],
        longitude,
        latitude,
    )

    # Temperature conversion
    if temperature_type == "potential":
        CT = gsw.CT_from_pt(SA, temp)

    elif temperature_type == "insitu":
        CT = gsw.CT_from_t(SA, temp, p[np.newaxis, :])

    else:
        raise ValueError("Unknown temperature_type")

    # Calculate seawater density
    rho = gsw.rho(SA, CT, p[np.newaxis, :])

    # Density anomalies relative to temporal mean
    rho_mean = np.nanmean(rho, axis=0)

    drho = rho - rho_mean[np.newaxis, :]

    # Integrate density anomalies from surface to sensor
    z_int = np.concatenate(
        [
            [0.0],
            z[(z > 0) & (z < depth_sensor)],
            [depth_sensor],
        ]
    )

    pressure_density = []

    for i in range(len(point.time)):

        valid = np.isfinite(drho[i, :])

        if valid.sum() < 2:
            pressure_density.append(np.nan)
            continue

        profile = np.interp(
            z_int,
            z[valid],
            drho[i, valid],
        )

        pressure_density.append(g * np.trapezoid(profile, z_int))

    pressure_density = np.asarray(pressure_density)

    # Sea surface height anomaly
    ssh_anomaly = ssh - np.nanmean(ssh)

    # Pressure contribution from SSH
    pressure_ssh = rho0 * g * ssh_anomaly

    # Total non-tidal ocean pressure anomaly
    dp_ocean = pressure_density + pressure_ssh

    result = xr.Dataset(
        {
            "dp_ocean": ("time", dp_ocean, {"units": "Pa"}),
            "dp_density": ("time", pressure_density, {"units": "Pa"}),
            "dp_ssh": ("time", pressure_ssh, {"units": "Pa"}),
            "ssh_anomaly": ("time", ssh_anomaly, {"units": "m"}),
            "equivalent_height": ("time", dp_ocean / (rho0 * g), {"units": "m"}),
        },
        coords={"time": point.time.values},
    )

    return result
