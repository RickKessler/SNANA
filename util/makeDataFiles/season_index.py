"""
season_index.py

Assigns an integer observing-season index to (RA, Dec, MJD) observations
at a fixed observing site, based on when each field is up at night.

Authors: Timothee Collinet-Adler, Richard Kessler
For queries email timotheea@uchicago.edu
"""

import numpy as np
from astropy.coordinates import EarthLocation, get_sun
from astropy.time import Time

SENTINEL = -9
DAYS_PER_YEAR = 365.25

SITES = {
    "rubin": EarthLocation.from_geodetic(
        lon=-70.7494, lat=-30.2446, height=2647.0),
    "ctio": EarthLocation.from_geodetic(
        lon=-70.8150, lat=-30.1652, height=2215.0),
    "mauna_kea": EarthLocation.from_geodetic(
        lon=-155.4681, lat=19.8283, height=4160.0),
}

# Alias, pointing at the same objects.
SITES["lsst"] = SITES["rubin"]


def get_season(data_dict, location, mjd_start, max_airmass=2.5, sun_alt_max=-12.0):
    """
    Assign an integer season index to each observation.

    Numbering is anchored at RA = 0: the season containing `mjd_start` at
    RA = 0 is season 1. Seasons count up from there and may be negative
    for observations before that band. Because the bands are diagonal, a
    single MJD can fall in different seasons at different RA.

    Parameters
    ----------
    data_dict : dict of {'RA': list, 'DEC': list, 'MJD': list}
        One entry per observation, in matching order across the three
        equal-length lists. RA and Dec in degrees; MJD as a timestamp
        including fractional day, not a night index.
    location : str
        Observing site. Accepted names are 'rubin' (or 'lsst'), 'ctio',
        and 'mauna_kea'.
    mjd_start : float
        Reference epoch defining season 1, per the numbering rule above.
    max_airmass : float, optional
        Airmass limit for the observability check, converted internally to
        a minimum altitude. Default 2.5 (altitude 23.6 deg), chosen loose
        so that marginal but genuine observations are not flagged.
    sun_alt_max : float, optional
        Solar altitude in degrees below which an observation is considered
        to be at night. Default -12.0 (nautical twilight).

    Returns
    -------
    numpy.ndarray of int, shape (N,)
        Season index per input observation, in input order. Entries failing
        the observability check are set to -9.

    Raises
    ------
    ValueError
        If `data_dict` is not a dict with the 'RA'/'DEC'/'MJD' keys, or
        their lists are not equal-length and non-empty; if RA, Dec, MJD,
        `mjd_start`, `max_airmass` or `sun_alt_max` fall outside plausible
        ranges; if `location` is an unknown site name.
    """

    # 1. Resolve the site.
    try:
        site = SITES[location.lower()]
    except (KeyError, AttributeError):
        raise ValueError(
            f"Unknown site {location!r}. Known sites: {sorted(SITES)}"
        ) from None

    # 2. Unpack the observation dict.
    if not isinstance(data_dict, dict):
        raise ValueError(f"data_dict must be a dict, got {type(data_dict).__name__}")
    required_keys = ("RA", "DEC", "MJD")
    missing = [k for k in required_keys if k not in data_dict]
    if missing:
        raise ValueError(f"data_dict is missing required keys: {missing}")
    ra = np.asarray(data_dict["RA"], dtype=float)
    dec = np.asarray(data_dict["DEC"], dtype=float)
    mjd = np.asarray(data_dict["MJD"], dtype=float)
    if ra.ndim != 1 or not (ra.shape == dec.shape == mjd.shape):
        raise ValueError(
            "RA, DEC, and MJD must be equal-length 1-D sequences, got "
            f"shapes {ra.shape}, {dec.shape}, {mjd.shape}"
        )
    if ra.size == 0:
        raise ValueError("data_dict must contain at least one observation")

    # 3. Range-check the observations.
    if np.any((ra < 0.0) | (ra >= 360.0)):
        raise ValueError("RA must be in [0, 360) degrees")
    if np.any(np.abs(dec) > 90.0):
        raise ValueError("Dec must be in [-90, 90] degrees")
    if not np.all(np.isfinite(mjd)):
        raise ValueError("MJD values must be finite")

    # 4. Range-check the scalar arguments.
    if not np.isfinite(mjd_start):
        raise ValueError(f"mjd_start must be finite, got {mjd_start!r}")
    if not (np.isfinite(max_airmass) and max_airmass > 1.0):
        raise ValueError(f"max_airmass must be a finite number > 1, got {max_airmass!r}")
    if not (np.isfinite(sun_alt_max) and -90.0 <= sun_alt_max <= 90.0):
        raise ValueError(f"sun_alt_max must be in [-90, 90] degrees, got {sun_alt_max!r}")

    # 5. Airmass limit as a minimum altitude.
    alt_min = np.degrees(np.arcsin(1.0 / max_airmass))

    # 6. Anchor the season count at the RA = 0 conjunction before mjd_start.
    t0 = _find_t0(mjd_start)

    # 7. Observability check.
    observable = _is_observable(mjd, ra, dec, site, alt_min, sun_alt_max)

    # 8. Season index. Shear the anchor from RA = 0 to this field's own
    # conjunction ladder, then count how many rungs have passed.
    season = np.floor((mjd - t0 - (ra / 360.0) * DAYS_PER_YEAR) / DAYS_PER_YEAR) + 1

    # 9. Sentinel where the observability check failed. Indices are computed
    #    for every row because masking after is cheaper than branching per row.
    return np.where(observable, season, SENTINEL).astype(int)


def _is_observable(mjd, ra, dec, site, alt_min, sun_alt_max):
    """
    Return a boolean mask marking observations as physically possible.

    An observation passes if the field is above `alt_min` and the Sun is
    below `sun_alt_max` at that instant. Both are evaluated at the exact
    time given by `mjd`, so the MJD must carry a fractional day.

    Altitude and solar position are computed with plain NumPy -- a
    low-precision solar ephemeris (Meeus/USNO, good to ~0.01 deg), a
    standard precession correction taking the input RA/Dec from J2000 to
    the equinox of date, and the standard sidereal-time / spherical-
    astronomy altitude formula -- instead of astropy's SkyCoord/AltAz
    transform pipeline, which was the runtime bottleneck at scale. This
    matches astropy's altitudes to within ~0.01 deg (nutation and
    aberration, the remaining uncorrected terms, are themselves at that
    level), well inside the margins `alt_min`/`sun_alt_max` are set to.

    Parameters
    ----------
    mjd, ra, dec : array_like of float
        Observation times (MJD) and field coordinates (degrees), equal length.
    site : astropy.coordinates.EarthLocation
        Already-resolved observing site.
    alt_min : float
        Minimum field altitude in degrees.
    sun_alt_max : float
        Maximum solar altitude in degrees.

    Returns
    -------
    numpy.ndarray of bool
    """

    jd = np.asarray(mjd, dtype=float) + 2400000.5
    n = jd - 2451545.0  # days since J2000.0
    t = n / 36525.0  # Julian centuries since J2000.0
    lat = np.radians(site.lat.deg)
    lon_deg = site.lon.deg

    # Low-precision solar position (apparent geocentric RA/Dec, of date).
    mean_lon = np.radians((280.460 + 0.9856474 * n) % 360.0)
    mean_anom = np.radians((357.528 + 0.9856003 * n) % 360.0)
    eclip_lon = (
        mean_lon
        + np.radians(1.915) * np.sin(mean_anom)
        + np.radians(0.020) * np.sin(2.0 * mean_anom)
    )
    obliquity = np.radians(23.439 - 0.0000004 * n)

    sun_ra = np.arctan2(np.cos(obliquity) * np.sin(eclip_lon), np.cos(eclip_lon))
    sun_dec = np.arcsin(np.sin(obliquity) * np.sin(eclip_lon))

    # Precess the field from J2000 (the input RA/Dec's frame) to the
    # equinox of date, since sidereal time is measured in that frame.
    arcsec = np.pi / (180.0 * 3600.0)
    zeta = (2306.2181 * t + 0.30188 * t**2 + 0.017998 * t**3) * arcsec
    z = (2306.2181 * t + 1.09468 * t**2 + 0.018203 * t**3) * arcsec
    theta = (2004.3109 * t - 0.42665 * t**2 - 0.041833 * t**3) * arcsec

    ra0 = np.radians(ra)
    dec0 = np.radians(dec)
    a = np.cos(dec0) * np.sin(ra0 + zeta)
    b = np.cos(theta) * np.cos(dec0) * np.cos(ra0 + zeta) - np.sin(theta) * np.sin(dec0)
    c = np.sin(theta) * np.cos(dec0) * np.cos(ra0 + zeta) + np.cos(theta) * np.sin(dec0)
    field_ra = np.arctan2(a, b) + z
    field_dec = np.arcsin(c)

    # Greenwich mean sidereal time -> local sidereal time at the site.
    gmst_deg = (
        280.46061837
        + 360.98564736629 * n
        + 0.000387933 * t**2
        - t**3 / 38710000.0
    ) % 360.0
    lst = np.radians((gmst_deg + lon_deg) % 360.0)

    def altitude_deg(obj_ra_rad, obj_dec_rad):
        hour_angle = lst - obj_ra_rad
        sin_alt = (
            np.sin(obj_dec_rad) * np.sin(lat)
            + np.cos(obj_dec_rad) * np.cos(lat) * np.cos(hour_angle)
        )
        return np.degrees(np.arcsin(sin_alt))

    field_alt = altitude_deg(field_ra, field_dec)
    sun_alt = altitude_deg(sun_ra, sun_dec)

    return (sun_alt < sun_alt_max) & (field_alt > alt_min)


def _find_t0(mjd_start):
    """
    Return the most recent MJD at or before `mjd_start` when the Sun's RA
    was zero. This anchors the season count: the band beginning here is
    season 1 at RA = 0.

    The Sun's RA climbs from 0 to 360 degrees once per year and wraps back
    to zero. The wrap is located on a daily grid, then
    refined by linear interpolation across the day it falls in.
    """
    grid = np.arange(mjd_start - 370.0, mjd_start + 1.0, 1.0)
    ra = get_sun(Time(grid, format="mjd")).ra.deg

    # The only place the Sun's RA decreases is the 360 -> 0 wrap.
    drops = np.flatnonzero(np.diff(ra) < 0)
    i = drops[-1]

    # Interpolate across the wrap: RA continues past 360 into the next day.
    frac = (360.0 - ra[i]) / (ra[i + 1] + 360.0 - ra[i])
    return grid[i] + frac
