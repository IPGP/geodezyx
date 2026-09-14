"""
Reading routine for RAF20.tac grid file (French geoid model NGF-IGN69 in RGF93).

Based on the RAF20.tac format specification:
- Grid extent: 42° ≤ φ ≤ 51.5°, -5.5° ≤ λ ≤ 8.5°
- Grid step: 0.025° in latitude, 0.0333° in longitude
- Direct readable format by Circé 5.2 software
- Data organized in blocks: each block has 381 latitude lines with 10 longitude values each
"""

import numpy as np
from typing import Dict, Tuple, Any


def read_raf20_tac(filepath: str) -> Dict[str, Any]:
    """
    Read RAF20.tac geoid grid file (French geoid model).

    File Structure:
    - Header line with grid parameters
    - Data organized in blocks: each block has ~381 lines (one per latitude)
    - Each line: 10 geoid height values (+ quality codes) for different longitudes
    - Total: 42 blocks × 10 longitudes + partial block = 421 longitude indices

    Parameters
    ----------
    filepath : str
        Path to the RAF20.tac file

    Returns
    -------
    dict
        Dictionary containing:
        - Grid metadata: lon_min, lon_max, lat_min, lat_max, dlon, dlat, description
        - Grid dimensions: n_lon_max, n_lat_max, n_actual_points
        - Data: points (list of tuples), point_dict (lookup dict)
        - Coordinates: grid_lons, grid_lats (sorted unique values)
    """

    result = {}

    with open(filepath, 'r', encoding='utf-8') as f:
        # Parse header line
        header_line = f.readline().strip()
        header_parts = header_line.split()

        # Extract grid parameters
        lon_min = float(header_parts[0])
        lon_max = float(header_parts[1])
        lat_min = float(header_parts[2])
        lat_max = float(header_parts[3])
        dlon = float(header_parts[4])
        dlat = float(header_parts[5])

        # Extract description (skip numeric fields 0-8, join rest)
        description = ' '.join(header_parts[9:])

        # Calculate expected grid dimensions
        n_lon_max = int(np.round((lon_max - lon_min) / dlon)) + 1
        n_lat_max = int(np.round((lat_max - lat_min) / dlat)) + 1

        # Store metadata
        result.update({
            'header': header_line,
            'lon_min': lon_min, 'lon_max': lon_max,
            'lat_min': lat_min, 'lat_max': lat_max,
            'dlon': dlon, 'dlat': dlat,
            'description': description,
            'n_lon_max': n_lon_max,
            'n_lat_max': n_lat_max,
        })

        # Read data organized in longitude blocks
        points = []
        point_dict = {}
        unique_lons = set()
        unique_lats = set()

        lon_idx = 0  # Current 10-longitude block index
        lat_idx = 0  # Current latitude within block

        for line in f:
            if not line.strip():
                continue

            tokens = line.split()
            n_pairs = len(tokens) // 2

            if n_pairs > 0:
                # Calculate current latitude
                lat = lat_min + lat_idx * dlat
                unique_lats.add(lat)

                # Parse geoid height values (each pair is: value, quality_code)
                for i in range(n_pairs):
                    try:
                        value = float(tokens[2*i])
                        qual = int(tokens[2*i + 1])

                        # Calculate longitude for this sub-value
                        # Each line contains 10 sub-longitudes within the block
                        lon = lon_min + (lon_idx * 10 + i) * dlon
                        unique_lons.add(lon)

                        # Store point
                        points.append((lon, lat, value, qual))

                        # Add to dictionary with rounded coordinates as key
                        key = (round(lon, 10), round(lat, 10))
                        point_dict[key] = (value, qual)

                    except (ValueError, IndexError):
                        break

                lat_idx += 1

                # Move to next longitude block after reading all latitudes
                if lat_idx >= n_lat_max:
                    lat_idx = 0
                    lon_idx += 1

        # Store data
        result.update({
            'points': points,
            'point_dict': point_dict,
            'grid_lons': np.array(sorted(unique_lons)),
            'grid_lats': np.array(sorted(unique_lats)),
            'n_actual_points': len(points),
        })

    return result


def get_raf20_value_at_point(raf20_data: Dict[str, Any], lon: float, lat: float) -> Tuple[float, int, bool]:
    """
    Get geoid height value at a specific geographic point.

    Tries to find an exact matching grid point first, then uses nearest neighbor.

    Parameters
    ----------
    raf20_data : dict
        Dictionary returned by read_raf20_tac()
    lon : float
        Longitude in degrees
    lat : float
        Latitude in degrees

    Returns
    -------
    tuple
        (geoid_height, quality_code, valid_flag)
        - geoid_height: computed geoid height in meters
        - quality_code: quality indicator from grid
        - valid_flag: True if a value was found
    """

    dlon = raf20_data['dlon']
    dlat = raf20_data['dlat']
    point_dict = raf20_data['point_dict']

    # Try to find exact grid point first
    lon_rounded = round(lon / dlon) * dlon
    lat_rounded = round(lat / dlat) * dlat
    key = (round(lon_rounded, 10), round(lat_rounded, 10))

    if key in point_dict:
        value, quality = point_dict[key]
        return value, quality, True

    # Use nearest neighbor search if exact point not found
    min_dist = float('inf')
    nearest_value = np.nan
    nearest_quality = -1

    for (grid_lon, grid_lat), (value, quality) in point_dict.items():
        # Euclidean distance in degree space
        dist = (grid_lon - lon)**2 + (grid_lat - lat)**2
        if dist < min_dist:
            min_dist = dist
            nearest_value = value
            nearest_quality = quality

    return (nearest_value, nearest_quality, True) if min_dist < float('inf') else (np.nan, -1, False)


def print_raf20_info(raf20_data: Dict[str, Any]) -> None:
    """Print summary information about the RAF20.tac grid."""

    print("RAF20 Geoid Grid File Information")
    print("=" * 70)
    print(f"Description: {raf20_data['description']}")
    print(f"\nGrid Coverage:")
    print(f"  Latitude:  {raf20_data['lat_min']:7.2f}° to {raf20_data['lat_max']:7.2f}°")
    print(f"  Longitude: {raf20_data['lon_min']:7.2f}° to {raf20_data['lon_max']:7.2f}°")
    print(f"\nGrid Spacing:")
    print(f"  Latitude step:  {raf20_data['dlat']:.6f}°")
    print(f"  Longitude step: {raf20_data['dlon']:.6f}°")
    print(f"\nGrid Dimensions:")
    print(f"  Max grid size: {raf20_data['n_lon_max']} longitude × {raf20_data['n_lat_max']} latitude points")
    print(f"  Actual points: {raf20_data['n_actual_points']} ({100.0*raf20_data['n_actual_points']/(raf20_data['n_lon_max']*raf20_data['n_lat_max']):.1f}% coverage)")
    print(f"  Unique longitudes: {len(raf20_data['grid_lons'])}")
    print(f"  Unique latitudes: {len(raf20_data['grid_lats'])}")

    values = np.array([v for _, _, v, _ in raf20_data['points']])
    print(f"\nGeoid Height Statistics:")
    print(f"  Minimum: {np.min(values):8.4f} m")
    print(f"  Maximum: {np.max(values):8.4f} m")
    print(f"  Mean:    {np.mean(values):8.4f} m")
    print(f"  Std dev: {np.std(values):8.4f} m")


if __name__ == "__main__":
    # Example usage
    filepath = "/home/sakic/Downloads/RAF20.tac"

    print("Reading RAF20.tac file...\n")
    raf20 = read_raf20_tac(filepath)

    print_raf20_info(raf20)

    # Example: Get value at Paris
    print("\n" + "=" * 70)
    print("Example: Geoid height at Paris (2.35°E, 48.85°N)")
    value, quality, valid = get_raf20_value_at_point(raf20, 2.35, 48.85)
    if valid:
        print(f"  Geoid height: {value:.4f} m")
        print(f"  Quality code: {quality}")
    else:
        print("  Point is outside grid bounds")

    # Display first few points
    print("\n" + "=" * 70)
    print("First 10 points in the grid:")
    for i, (lon, lat, value, quality) in enumerate(raf20['points'][:10]):
        print(f"  {i+1:2d}. ({lat:8.4f}°, {lon:8.4f}°): {value:8.4f} m (quality: {quality})")

