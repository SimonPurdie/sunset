import numpy as np
from skyfield.api import load
from skyfield.toposlib import Topos
from skyfield.errors import EphemerisRangeError
from typing import Optional
import warnings

from sunset.models import Body, Observer, SOLAR_SYSTEM_BODIES


def _parse_iso_utc(utc_time: str):
    """
    Parse ISO-8601 UTC timestamp and return components for skyfield.

    Args:
        utc_time: ISO-8601 UTC timestamp (e.g., "2026-01-20T12:00:00Z")

    Returns:
        Tuple of (year, month, day, hour, minute, second)
    """
    from datetime import datetime

    if utc_time.endswith("Z"):
        utc_time = utc_time[:-1] + "+00:00"

    dt = datetime.fromisoformat(utc_time)
    return (dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)


def resolve_sunset_location(
    utc_time: str,
    body_id: str,
    random_seed: Optional[int] = None,
) -> Observer:
    """
    Find an observer location on the given body where the sun is at sunset
    (solar elevation between -1.5° and +0.5°).

    Args:
        utc_time: ISO-8601 UTC timestamp
        body_id: Solar System body identifier (e.g., "earth", "mars")
        random_seed: Optional seed for deterministic random selection

    Returns:
        Observer dataclass with location and solar geometry

    Raises:
        ValueError: If body_id is not recognized or no sunset found
    """
    if body_id not in SOLAR_SYSTEM_BODIES:
        raise ValueError(f"Unknown body: {body_id}")

    body = SOLAR_SYSTEM_BODIES[body_id]

    ephemeris = load("de421.bsp")
    sun = ephemeris["sun"]

    ts = load.timescale()
    t = ts.utc(*_parse_iso_utc(utc_time))

    latitude, longitude, solar_elevation, sun_direction, solar_angular_diameter = (
        _find_sunset_coordinates(t, ephemeris, sun, body_id, body.radius_m, random_seed)
    )

    return Observer(
        latitude=latitude,
        longitude=longitude,
        altitude_m=0.0,
        body=body,
        solar_elevation_deg=solar_elevation,
        sun_direction=sun_direction,
        solar_angular_diameter_deg=solar_angular_diameter,
    )


def _find_sunset_coordinates(
    t, ephemeris, sun, body_id: str, body_radius_m: float, random_seed: Optional[int]
) -> tuple[float, float, float, tuple[float, float, float], float]:
    """
    Find latitude and longitude where solar elevation is between -1.5° and +0.5°.

    Uses random sampling and selects a point deterministically based on seed.

    Returns:
        (latitude, longitude, solar_elevation, sun_direction, solar_angular_diameter)
    """
    rng = np.random.default_rng(random_seed)

    best_elevation = -99.0
    best_latitude = 0.0
    best_longitude = 0.0
    best_sun_direction = (0.0, 0.0, 0.0)
    best_solar_angular_diameter = 0.0

    planet = ephemeris[body_id]

    for _ in range(500):
        test_lat = rng.uniform(-90, 90)
        test_lon = rng.uniform(-180, 180)

        observer = planet + Topos(latitude_degrees=test_lat, longitude_degrees=test_lon)
        astrometric = sun.at(t).observe(observer)

        alt, az, distance = _get_alt_az(astrometric, test_lat, test_lon)
        elevation = float(alt.degrees)

        if -1.5 <= elevation <= 0.5:
            if abs(elevation) < abs(best_elevation):
                best_elevation = elevation
                best_latitude = test_lat
                best_longitude = test_lon
                best_sun_direction = _az_alt_to_vector(
                    float(az.degrees), float(alt.degrees)
                )
                best_solar_angular_diameter = _calculate_solar_angular_diameter(
                    float(distance.km)
                )

        if abs(elevation) < 0.01:
            break

    if abs(best_elevation) > 2.0:
        raise ValueError(
            f"Could not find sunset location for {body_id} at {t.utc_strftime()}"
        )

    return (
        best_latitude,
        best_longitude,
        best_elevation,
        best_sun_direction,
        best_solar_angular_diameter,
    )


def _get_alt_az(astrometric, lat_deg: float, lon_deg: float):
    """
    Get altitude and azimuth from astrometric position, handling skyfield API limitations.

    This works around a bug in skyfield's relativistic calculations that can cause
    .apparent() to fail with EphemerisRangeError.
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            apparent = astrometric.apparent()
        return apparent.altaz()
    except EphemerisRangeError:
        return _alt_az_from_radec(astrometric, lat_deg, lon_deg)


def _alt_az_from_radec(astrometric, lat_deg: float, lon_deg: float):
    """
    Calculate altitude and azimuth from RA/Dec coordinates.

    This is a fallback when .apparent() fails due to ephemeris range issues.
    """
    from skyfield.units import Angle
    import math

    ra, dec, distance = astrometric.radec(epoch="date")

    sidereal_time = _get_local_sidereal_time(astrometric.t.utc_strftime(), lon_deg)
    hour_angle = sidereal_time - float(ra.hours)

    lat_rad = math.radians(lat_deg)
    ha_rad = math.radians(hour_angle * 15.0)
    dec_rad = float(dec.radians)

    alt_rad = math.asin(
        math.sin(lat_rad) * math.sin(dec_rad)
        + math.cos(lat_rad) * math.cos(dec_rad) * math.cos(ha_rad)
    )
    az_rad = math.atan2(
        -math.sin(ha_rad),
        math.cos(lat_rad) * math.tan(dec_rad) - math.sin(lat_rad) * math.cos(ha_rad),
    )

    altitude = Angle(degrees=math.degrees(alt_rad))
    azimuth = Angle(degrees=(math.degrees(az_rad) % 360.0))

    return (altitude, azimuth, distance)


def _get_local_sidereal_time(utc_str, lon_deg: float) -> float:
    """
    Calculate Local Sidereal Time in hours.

    Simple approximation - for more accurate results would use more precise formulas.
    """
    from datetime import datetime
    import re

    match = re.match(
        r"(\d{4})-(\d{2})-(\d{2})\s+(\d{2}):(\d{2}):(\d{2})\s+UTC", utc_str
    )
    if match:
        year, month, day, hour, minute, second = map(int, match.groups())
    else:
        raise ValueError(f"Could not parse UTC string: {utc_str}")

    jd0 = _julian_day(year, month, day)
    jd = jd0 + (hour + minute / 60.0 + second / 3600.0) / 24.0

    t = (jd - 2451545.0) / 36525.0
    gmst = (
        280.46061837
        + 360.98564736629 * (jd - 2451545.0)
        + 0.000387933 * t * t
        - t * t * t / 38710000.0
    ) % 360.0

    gmst_hours = gmst / 15.0
    lst_hours = (gmst_hours + lon_deg / 15.0) % 24.0

    return lst_hours


def _julian_day(year: int, month: int, day: int) -> float:
    """Calculate Julian Day number."""
    if month <= 2:
        year -= 1
        month += 12

    a = year // 100
    b = 2 - a + a // 4

    return int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + b - 1524.5


def _calculate_solar_angular_diameter(distance_km: float) -> float:
    """Calculate solar angular diameter in degrees based on distance."""
    solar_radius_km = 696_340
    angular_diameter_rad = 2 * np.arctan(solar_radius_km / distance_km)
    return np.degrees(angular_diameter_rad)


def _az_alt_to_vector(az_deg: float, alt_deg: float) -> tuple[float, float, float]:
    """
    Convert azimuth/altitude to unit direction vector in ENU (East-North-Up) coordinates.

    Args:
        az_deg: Azimuth angle in degrees (clockwise from north)
        alt_deg: Altitude angle in degrees (elevation above horizon)

    Returns:
        Unit vector as (east, north, up) tuple
    """
    alt_rad = np.radians(alt_deg)
    az_rad = np.radians(az_deg)

    east = np.cos(alt_rad) * np.sin(az_rad)
    north = np.cos(alt_rad) * np.cos(az_rad)
    up = np.sin(alt_rad)

    return (east, north, up)
