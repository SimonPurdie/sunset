import pytest
from sunset.models import (
    Body,
    Observer,
    AtmosphericProfile,
    NullAtmosphere,
    SpectralRadiance,
    Scene,
    SOLAR_SYSTEM_BODIES,
    get_atmosphere,
)
import numpy as np


def test_body_dataclass():
    body = SOLAR_SYSTEM_BODIES["earth"]
    assert body.name == "Earth"
    assert body.radius_m == 6_371_000
    assert body.has_atmosphere is True


def test_airless_bodies():
    mercury = SOLAR_SYSTEM_BODIES["mercury"]
    moon = SOLAR_SYSTEM_BODIES["moon"]
    assert mercury.has_atmosphere is False
    assert moon.has_atmosphere is False


def test_all_bodies_have_required_fields():
    for body_id, body in SOLAR_SYSTEM_BODIES.items():
        assert isinstance(body.name, str)
        assert body.radius_m > 0
        assert isinstance(body.has_atmosphere, bool)
        assert len(body.ground_color_base) == 3
        assert 0.0 <= body.ground_color_variation <= 1.0


def test_observer_creation():
    body = SOLAR_SYSTEM_BODIES["earth"]
    observer = Observer(
        latitude=40.0,
        longitude=-74.0,
        altitude_m=100.0,
        body=body,
        solar_elevation_deg=0.0,
        sun_direction=(0.0, 1.0, 0.0),
        solar_angular_diameter_deg=0.5,
    )
    assert observer.latitude == 40.0
    assert observer.longitude == -74.0
    assert observer.altitude_m == 100.0
    assert observer.body == body
    assert observer.solar_elevation_deg == 0.0
    assert observer.sun_direction == (0.0, 1.0, 0.0)
    assert observer.solar_angular_diameter_deg == 0.5


def test_null_atmosphere_for_airless_bodies():
    mercury_atmosphere = get_atmosphere("mercury", altitude_m=0)
    assert isinstance(mercury_atmosphere, NullAtmosphere)
    assert mercury_atmosphere.body_name == "Mercury"

    moon_atmosphere = get_atmosphere("moon", altitude_m=0)
    assert isinstance(moon_atmosphere, NullAtmosphere)
    assert moon_atmosphere.body_name == "Moon"


def test_atmosphere_implemented_for_planets():
    earth_atmosphere = get_atmosphere("earth", altitude_m=0)
    assert isinstance(earth_atmosphere, AtmosphericProfile)
    assert earth_atmosphere.body_name == "Earth"

    mars_atmosphere = get_atmosphere("mars", altitude_m=0)
    assert isinstance(mars_atmosphere, AtmosphericProfile)
    assert mars_atmosphere.body_name == "Mars"

    venus_atmosphere = get_atmosphere("venus", altitude_m=0)
    assert isinstance(venus_atmosphere, AtmosphericProfile)
    assert venus_atmosphere.body_name == "Venus"

    titan_atmosphere = get_atmosphere("titan", altitude_m=0)
    assert isinstance(titan_atmosphere, AtmosphericProfile)
    assert titan_atmosphere.body_name == "Titan"


def test_spectral_radiance_validation():
    wavelengths = np.array([400, 500, 600, 700])
    radiance = np.array([1.0, 2.0, 3.0, 4.0])
    skylight = np.array([0.5, 1.0, 1.5, 2.0])

    spectral = SpectralRadiance(
        wavelengths=wavelengths,
        radiance=radiance,
        skylight_contribution=skylight,
    )
    assert len(spectral.wavelengths) == 4
    assert len(spectral.radiance) == 4
    assert len(spectral.skylight_contribution) == 4


def test_spectral_radiance_mismatched_lengths():
    wavelengths = np.array([400, 500, 600, 700])
    radiance = np.array([1.0, 2.0, 3.0])
    skylight = np.array([0.5, 1.0, 1.5, 2.0])

    with pytest.raises(
        ValueError, match="wavelengths and radiance must have same length"
    ):
        SpectralRadiance(
            wavelengths=wavelengths,
            radiance=radiance,
            skylight_contribution=skylight,
        )


def test_scene_creation():
    body = SOLAR_SYSTEM_BODIES["earth"]
    scene = Scene(
        body=body,
        utc_time="2026-01-20T12:00:00Z",
        latitude=40.0,
        longitude=-74.0,
        altitude_m=0.0,
        random_seed=12345,
        renderer_id="0.1.0",
    )
    assert scene.body == body
    assert scene.utc_time == "2026-01-20T12:00:00Z"
    assert scene.latitude == 40.0
    assert scene.longitude == -74.0
    assert scene.altitude_m == 0.0
    assert scene.random_seed == 12345
    assert scene.renderer_id == "0.1.0"
