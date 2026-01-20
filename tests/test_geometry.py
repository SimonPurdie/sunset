import pytest
import numpy as np
from sunset.geometry.resolver import resolve_sunset_location
from sunset.models import SOLAR_SYSTEM_BODIES


def test_solar_elevation_bounds():
    """Test that returned solar elevation is within [-1.5°, +0.5°]."""
    utc_time = "2026-01-20T18:00:00Z"
    body_id = "earth"
    random_seed = 12345

    observer = resolve_sunset_location(utc_time, body_id, random_seed)

    assert -1.5 <= observer.solar_elevation_deg <= 0.5


def test_solar_elevation_bounds_mars():
    """Test solar elevation bounds for Mars."""
    pytest.skip(
        "Skyfield's Topos is Earth-specific. Non-Earth planets require PlanetaryConstants "
        "with binary rotation data files (.bcp) that need to be downloaded separately. "
        "Implementation in progress - see docs/BREADCRUMBS.md for details."
    )
    utc_time = "2026-01-20T12:00:00Z"
    body_id = "mars"
    random_seed = 54321

    observer = resolve_sunset_location(utc_time, body_id, random_seed)

    assert -1.5 <= observer.solar_elevation_deg <= 0.5


def test_coordinates_valid():
    """Test that coordinates are valid for the body."""
    utc_time = "2026-01-20T18:00:00Z"
    body_id = "earth"
    random_seed = 98765

    observer = resolve_sunset_location(utc_time, body_id, random_seed)

    assert -90 <= observer.latitude <= 90
    assert -180 <= observer.longitude <= 180


def test_determinism():
    """Test that same seed produces same result."""
    utc_time = "2026-01-20T18:00:00Z"
    body_id = "earth"
    random_seed = 11111

    observer1 = resolve_sunset_location(utc_time, body_id, random_seed)
    observer2 = resolve_sunset_location(utc_time, body_id, random_seed)

    assert observer1.latitude == pytest.approx(observer2.latitude, abs=1e-6)
    assert observer1.longitude == pytest.approx(observer2.longitude, abs=1e-6)
    assert observer1.solar_elevation_deg == pytest.approx(
        observer2.solar_elevation_deg, abs=1e-6
    )
    assert observer1.sun_direction == observer2.sun_direction


def test_determinism_different_seeds():
    """Test that different seeds produce different results."""
    utc_time = "2026-01-20T18:00:00Z"
    body_id = "earth"

    observer1 = resolve_sunset_location(utc_time, body_id, 22222)
    observer2 = resolve_sunset_location(utc_time, body_id, 33333)

    latitude_diff = abs(observer1.latitude - observer2.latitude)
    longitude_diff = abs(observer1.longitude - observer2.longitude)

    assert latitude_diff > 0.001 or longitude_diff > 0.001


def test_sun_direction_is_unit_vector():
    """Test that sun direction is a normalized vector."""
    utc_time = "2026-01-20T18:00:00Z"
    body_id = "earth"
    random_seed = 44444

    observer = resolve_sunset_location(utc_time, body_id, random_seed)

    east, north, up = observer.sun_direction
    magnitude = np.sqrt(east**2 + north**2 + up**2)

    assert magnitude == pytest.approx(1.0, abs=1e-6)


def test_solar_angular_diameter_positive():
    """Test that solar angular diameter is positive and reasonable."""
    utc_time = "2026-01-20T18:00:00Z"
    body_id = "earth"
    random_seed = 55555

    observer = resolve_sunset_location(utc_time, body_id, random_seed)

    assert 0.1 < observer.solar_angular_diameter_deg < 2.0


def test_unknown_body_raises_error():
    """Test that unknown body raises ValueError."""
    utc_time = "2026-01-20T18:00:00Z"
    body_id = "unknown_planet"
    random_seed = 66666

    with pytest.raises(ValueError, match="Unknown body"):
        resolve_sunset_location(utc_time, body_id, random_seed)


def test_body_association():
    """Test that observer is correctly associated with body."""
    utc_time = "2026-01-20T18:00:00Z"
    body_id = "earth"
    random_seed = 77777

    observer = resolve_sunset_location(utc_time, body_id, random_seed)

    assert observer.body.name == "Earth"
    assert observer.body == SOLAR_SYSTEM_BODIES["earth"]


def test_mars_sunset():
    """Test that Mars sunset can be found."""
    pytest.skip(
        "Skyfield's Topos is Earth-specific. Non-Earth planets require PlanetaryConstants "
        "with binary rotation data files (.bcp) that need to be downloaded separately. "
        "Implementation in progress - see docs/BREADCRUMBS.md for details."
    )
    utc_time = "2026-01-20T12:00:00Z"
    body_id = "mars"
    random_seed = 88888

    observer = resolve_sunset_location(utc_time, body_id, random_seed)

    assert -1.5 <= observer.solar_elevation_deg <= 0.5
    assert observer.body.name == "Mars"


def test_venus_sunset():
    """Test that Venus sunset can be found."""
    pytest.skip(
        "Skyfield's Topos is Earth-specific. Non-Earth planets require PlanetaryConstants "
        "with binary rotation data files (.bcp) that need to be downloaded separately. "
        "Implementation in progress - see docs/BREADCRUMBS.md for details."
    )
    utc_time = "2026-01-20T00:00:00Z"
    body_id = "venus"
    random_seed = 99999

    observer = resolve_sunset_location(utc_time, body_id, random_seed)

    assert -1.5 <= observer.solar_elevation_deg <= 0.5
    assert observer.body.name == "Venus"


def test_mercury_sunset():
    """Test that Mercury sunset can be found."""
    pytest.skip(
        "Skyfield's Topos is Earth-specific. Non-Earth planets require PlanetaryConstants "
        "with binary rotation data files (.bcp) that need to be downloaded separately. "
        "Implementation in progress - see docs/BREADCRUMBS.md for details."
    )
    utc_time = "2026-01-20T06:00:00Z"
    body_id = "mercury"
    random_seed = 10101

    observer = resolve_sunset_location(utc_time, body_id, random_seed)

    assert -1.5 <= observer.solar_elevation_deg <= 0.5
    assert observer.body.name == "Mercury"
