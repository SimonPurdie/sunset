import numpy as np
from skyfield.api import load, PlanetaryConstants
from skyfield.toposlib import Topos
from skyfield.errors import EphemerisRangeError
from typing import Optional
import warnings

from sunset.models import Body, Observer, SOLAR_SYSTEM_BODIES
from sunset.errors import BodyNotFoundError, SunsetNotFoundError, TimeParseError


def _parse_iso_utc(utc_time: str):
    """
    Parse ISO-8601 UTC timestamp and return components for skyfield.

    Args:
        utc_time: ISO-8601 UTC timestamp (e.g., "2026-01-20T12:00:00Z")

    Returns:
        Tuple of (year, month, day, hour, minute, second)

    Raises:
        TimeParseError: If utc_time cannot be parsed
    """
    from datetime import datetime
    from sunset.errors import TimeParseError

    if utc_time.endswith("Z"):
        utc_time = utc_time[:-1] + "+00:00"

    try:
        dt = datetime.fromisoformat(utc_time)
    except ValueError:
        raise TimeParseError(utc_time)

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
        raise BodyNotFoundError(body_id, list(SOLAR_SYSTEM_BODIES.keys()))

    body = SOLAR_SYSTEM_BODIES[body_id]

    is_titan = body_id == "titan"
    if is_titan:
        ephemeris = load("de440.bsp")
        saturn = ephemeris[6]
        sun = ephemeris["sun"]
    else:
        ephemeris = load("de421.bsp")
        sun = ephemeris["sun"]
        saturn = None

    ts = load.timescale()
    t = ts.utc(*_parse_iso_utc(utc_time))

    latitude, longitude, solar_elevation, sun_direction, solar_angular_diameter = (
        _find_sunset_coordinates(
            t, ephemeris, sun, body_id, body.radius_m, random_seed, saturn
        )
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


def _get_titan_position(t, saturn):
    """
    Compute Titan's position relative to Saturn using orbital elements.

    Titan's orbital elements (J2000):
    - Semi-major axis (a): 1,221,870 km
    - Eccentricity (e): 0.0288
    - Inclination (i): 0.34854°
    - Longitude of ascending node (Ω): 168.469°
    - Argument of periapsis (ω): 179.09°
    - Mean anomaly at epoch (M0): 173.795°
    - Orbital period: 15.945 days

    Args:
        t: Skyfield Time object
        saturn: Skyfield Saturn object

    Returns:
        Skyfield position object for Titan
    """
    a = 1221870.0  # Semi-major axis in km
    e = 0.0288  # Eccentricity
    i = np.radians(0.34854)  # Inclination in radians
    omega = np.radians(168.469)  # Longitude of ascending node in radians
    w = np.radians(179.09)  # Argument of periapsis in radians
    M0 = np.radians(173.795)  # Mean anomaly at epoch in radians
    n = 2 * np.pi / (15.945 * 86400.0)  # Mean motion in rad/sec

    # Calculate time since J2000 epoch in seconds
    t_j2000 = load.timescale().utc(2000, 1, 1, 12, 0, 0)
    t_days = t.tai - t_j2000.tai
    dt = float(t_days * 86400.0)  # Convert days to seconds

    # Calculate current mean anomaly
    M = M0 + n * dt

    # Solve Kepler's equation for eccentric anomaly E using Newton-Raphson
    E = M
    for _ in range(10):
        E = M + e * np.sin(E)

    # Calculate true anomaly
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))

    # Calculate radius
    r = a * (1 - e * np.cos(E))

    # Calculate position in orbital plane
    x_orb = r * np.cos(nu)
    y_orb = r * np.sin(nu)

    # Rotate to ICRS coordinates
    cos_O = np.cos(omega)
    sin_O = np.sin(omega)
    cos_w = np.cos(w)
    sin_w = np.sin(w)
    cos_i = np.cos(i)
    sin_i = np.sin(i)

    # 3D rotation matrices
    x_ecl = x_orb * (cos_O * cos_w - sin_O * sin_w * cos_i) - y_orb * (
        cos_O * sin_w + sin_O * cos_w * cos_i
    )
    y_ecl = x_orb * (sin_O * cos_w + cos_O * sin_w * cos_i) + y_orb * (
        cos_O * cos_w * cos_i - sin_O * sin_w
    )
    z_ecl = x_orb * (sin_w * sin_i) + y_orb * (cos_w * sin_i)

    # Get Saturn's position
    saturn_pos = saturn.at(t).position.km

    # Add Titan's offset to Saturn's position
    titan_pos = saturn_pos + np.array([x_ecl, y_ecl, z_ecl])

    # Create a mock position object that mimics Skyfield's position interface
    class TitanPosition:
        def __init__(self, pos_km):
            self._pos_km = np.array(pos_km, dtype=np.float64)

        def at(self, t):
            # Return a position object with .position attribute
            class Position:
                def __init__(self, pos_km):
                    self.position = type(
                        "obj",
                        (object,),
                        {
                            "km": np.array(pos_km, dtype=np.float64),
                            "au": pos_km / 149597870.7,
                        },
                    )

            return Position(self._pos_km)

    return TitanPosition(titan_pos)


def _create_planetary_observer(
    planet, lat_deg: float, lon_deg: float, pc, frame_name: str
):
    """
    Create a surface observer position on a planet using PlanetaryConstants.

    Args:
        planet: Skyfield planet object
        lat_deg: Latitude in degrees
        lon_deg: Longitude in degrees
        pc: PlanetaryConstants object with loaded data
        frame_name: Name of the reference frame (e.g., "MARS_IAU2000")

    Returns:
        Skyfield position object representing the surface location
    """
    frame = pc.build_frame_named(frame_name)
    return planet + pc.build_latlon_degrees(frame, lat_deg, lon_deg)


def _find_sunset_coordinates(
    t,
    ephemeris,
    sun,
    body_id: str,
    body_radius_m: float,
    random_seed: Optional[int],
    saturn=None,
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

    is_titan = body_id == "titan"
    if is_titan:
        planet = _get_titan_position(t, saturn)
    else:
        planet = ephemeris[body_id]

    is_earth = body_id == "earth"

    if is_earth:
        for _ in range(500):
            test_lat = rng.uniform(-90, 90)
            test_lon = rng.uniform(-180, 180)

            observer = planet + Topos(
                latitude_degrees=test_lat, longitude_degrees=test_lon
            )
            astrometric = sun.at(t).observe(observer)
            alt, az, distance = _get_alt_az(astrometric, test_lat, test_lon)
            elevation = float(alt.degrees)

            if -1.5 <= elevation <= 0.5:
                if abs(elevation) < abs(best_elevation):
                    best_elevation = elevation
                    best_latitude = test_lat
                    best_longitude = test_lon
                    best_sun_direction = _az_alt_to_vector(float(az.degrees), elevation)
                    best_solar_angular_diameter = _calculate_solar_angular_diameter(
                        float(distance.km)
                    )

            if abs(elevation) < 0.01:
                break
    else:
        (
            best_latitude,
            best_longitude,
            best_elevation,
            best_sun_direction,
            best_solar_angular_diameter,
        ) = _find_planetary_sunset_terminator(
            t, ephemeris, sun, planet, body_id, body_radius_m, random_seed
        )

    if abs(best_elevation) > 2.0:
        raise SunsetNotFoundError(body_id, t.utc_strftime())

    return (
        best_latitude,
        best_longitude,
        best_elevation,
        best_sun_direction,
        best_solar_angular_diameter,
    )


def _find_planetary_sunset_terminator(
    t,
    ephemeris,
    sun,
    planet,
    body_id: str,
    body_radius_m: float,
    random_seed: Optional[int],
) -> tuple[float, float, float, tuple[float, float, float], float]:
    """
    Find sunset location on a non-Earth planet using direct geometric approach.

    This method finds points on the terminator (circle where sun is on horizon)
    using pure geometry without requiring planet rotation data.

    Args:
        t: Skyfield Time object
        ephemeris: Skyfield ephemeris
        sun: Skyfield sun object
        planet: Skyfield planet object
        body_id: Planet identifier
        body_radius_m: Planet radius in meters
        random_seed: Random seed for deterministic selection

    Returns:
        (latitude, longitude, solar_elevation, sun_direction, solar_angular_diameter)
    """
    planet_position = planet.at(t)
    sun_vec = sun.at(t).position.km
    planet_center = planet_position.position.km

    rng = np.random.default_rng(random_seed)

    sun_from_planet = sun_vec - planet_center
    sun_direction_norm = sun_from_planet / np.linalg.norm(sun_from_planet)

    best_elevation = 99.0
    best_lat = 0.0
    best_lon = 0.0
    best_sun_direction = (0.0, 0.0, 0.0)
    best_solar_angular_diameter = 0.0

    for i in range(2000):
        theta = rng.uniform(0, 2 * np.pi)
        phi = rng.uniform(-np.pi / 2, np.pi / 2)

        radial_dir = np.array(
            [np.cos(phi) * np.cos(theta), np.cos(phi) * np.sin(theta), np.sin(phi)]
        )

        dot_with_sun = np.dot(radial_dir, sun_direction_norm)
        dot_with_sun = np.clip(dot_with_sun, -1, 1)
        elevation = 90.0 - np.degrees(np.arccos(dot_with_sun))

        if -1.5 <= elevation <= 0.5 and abs(elevation) < abs(best_elevation):
            best_elevation = elevation
            best_lat = np.degrees(phi)
            best_lon = (np.degrees(theta) + 180) % 360 - 180

            up = radial_dir
            north = np.array([0, 0, 1])

            if np.abs(np.dot(up, north)) > 0.99:
                east = np.array([1, 0, 0])
                north = np.array([0, 0, 1])
            else:
                east = np.cross(north, up)
                east = east / np.linalg.norm(east)
                north = np.cross(up, east)

            east_comp = float(np.dot(sun_direction_norm, east))
            north_comp = float(np.dot(sun_direction_norm, north))
            up_comp = float(np.dot(sun_direction_norm, up))
            best_sun_direction = (east_comp, north_comp, up_comp)

            sun_distance_km = np.linalg.norm(sun_from_planet)
            best_solar_angular_diameter = _calculate_solar_angular_diameter(
                sun_distance_km
            )

        if abs(elevation) < 0.01:
            break

    if abs(best_elevation) > 2.0:
        raise SunsetNotFoundError(body_id, t.utc_strftime())

    return (
        best_lat,
        best_lon,
        best_elevation,
        best_sun_direction,
        best_solar_angular_diameter,
    )


def _compute_enu_direction_from_icrs(
    to_sun_norm, radial, lat_deg: float, lon_deg: float
) -> tuple[float, float, float]:
    """
    Compute sun direction in ENU coordinates from ICRS vectors.

    Args:
        to_sun_norm: Unit vector from surface point to sun (ICRS)
        radial: Unit radial vector from planet center to surface point (ICRS)
        lat_deg: Latitude of the surface point
        lon_deg: Longitude of the surface point

    Returns:
        (east, north, up) unit vector
    """
    lat_rad = np.radians(lat_deg)
    lon_rad = np.radians(lon_deg)

    east_body = np.array([-np.sin(lon_rad), np.cos(lon_rad), 0.0])
    north_body = np.array(
        [
            -np.sin(lat_rad) * np.cos(lon_rad),
            -np.sin(lat_rad) * np.sin(lon_rad),
            np.cos(lat_rad),
        ]
    )
    up_body = np.array(
        [
            np.cos(lat_rad) * np.cos(lon_rad),
            np.cos(lat_rad) * np.sin(lon_rad),
            np.sin(lat_rad),
        ]
    )

    east_comp = float(np.dot(to_sun_norm, east_body))
    north_comp = float(np.dot(to_sun_norm, north_body))
    up_comp = float(np.dot(to_sun_norm, up_body))

    return (east_comp, north_comp, up_comp)


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
        raise TimeParseError(utc_str)

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
