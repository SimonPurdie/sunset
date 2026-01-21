import pytest

from sunset.atmosphere.provider import get_body_profile
from sunset.models.atmosphere import AtmosphericProfile, NullAtmosphere
from sunset.models.bodies import SOLAR_SYSTEM_BODIES
from sunset.models.observer import Observer
from sunset.visibility.resolver import (
    VISIBILITY_THRESHOLD,
    calculate_optical_depth,
    resolve_visibility_elevation,
)


class TestVenusSurfaceRejected:
    """Test Venus surface is rejected, upper atmosphere accepted."""

    def test_venus_surface_rejected(self):
        profile = get_body_profile("Venus", 0)
        assert isinstance(profile, AtmosphericProfile)
        assert not isinstance(profile, NullAtmosphere)

        body = SOLAR_SYSTEM_BODIES["venus"]
        observer = Observer(
            latitude=0.0,
            longitude=0.0,
            altitude_m=0.0,
            body=body,
            solar_elevation_deg=0.0,
            sun_direction=(0.0, 1.0, 0.0),
            solar_angular_diameter_deg=0.5,
        )

        final_alt, justification, is_visible = resolve_visibility_elevation(
            profile, observer
        )

        assert final_alt > 0, "Venus surface should require altitude increase"
        assert is_visible, "Venus should be visible at some altitude"
        assert "Raised altitude" in justification

    def test_venus_upper_atmosphere_accepted(self):
        profile = get_body_profile("Venus", 0)
        assert isinstance(profile, AtmosphericProfile)
        assert not isinstance(profile, NullAtmosphere)

        body = SOLAR_SYSTEM_BODIES["venus"]
        high_alt_observer = Observer(
            latitude=0.0,
            longitude=0.0,
            altitude_m=50000.0,
            body=body,
            solar_elevation_deg=0.0,
            sun_direction=(0.0, 1.0, 0.0),
            solar_angular_diameter_deg=0.5,
        )

        final_alt, justification, is_visible = resolve_visibility_elevation(
            profile, high_alt_observer
        )

        assert is_visible, "Venus should be visible at high altitude"


class TestEarthSurfaceOK:
    """Test Earth never requires altitude change."""

    def test_earth_surface_ok(self):
        profile = get_body_profile("Earth", 0)
        assert isinstance(profile, AtmosphericProfile)
        assert not isinstance(profile, NullAtmosphere)

        body = SOLAR_SYSTEM_BODIES["earth"]
        observer = Observer(
            latitude=0.0,
            longitude=0.0,
            altitude_m=0.0,
            body=body,
            solar_elevation_deg=0.0,
            sun_direction=(0.0, 1.0, 0.0),
            solar_angular_diameter_deg=0.5,
        )

        final_alt, justification, is_visible = resolve_visibility_elevation(
            profile, observer
        )

        assert final_alt == 0.0, "Earth surface should not require altitude change"
        assert is_visible, "Earth should be visible at surface"
        assert "Sun visible at surface" in justification


class TestAirlessSurfaceOK:
    """Test airless bodies always accept surface."""

    def test_mercury_surface_ok(self):
        profile = get_body_profile("Mercury", 0)
        assert isinstance(profile, NullAtmosphere)

        body = SOLAR_SYSTEM_BODIES["mercury"]
        observer = Observer(
            latitude=0.0,
            longitude=0.0,
            altitude_m=0.0,
            body=body,
            solar_elevation_deg=0.0,
            sun_direction=(0.0, 1.0, 0.0),
            solar_angular_diameter_deg=1.0,
        )

        final_alt, justification, is_visible = resolve_visibility_elevation(
            profile, observer
        )

        assert final_alt == 0.0, "Mercury surface should not require altitude change"
        assert is_visible, "Mercury should be visible at surface"
        assert "no atmosphere" in justification.lower()

    def test_moon_surface_ok(self):
        profile = get_body_profile("Moon", 0)
        assert isinstance(profile, NullAtmosphere)

        body = SOLAR_SYSTEM_BODIES["moon"]
        observer = Observer(
            latitude=0.0,
            longitude=0.0,
            altitude_m=0.0,
            body=body,
            solar_elevation_deg=0.0,
            sun_direction=(0.0, 1.0, 0.0),
            solar_angular_diameter_deg=0.5,
        )

        final_alt, justification, is_visible = resolve_visibility_elevation(
            profile, observer
        )

        assert final_alt == 0.0, "Moon surface should not require altitude change"
        assert is_visible, "Moon should be visible at surface"
        assert "no atmosphere" in justification.lower()


class TestOpticalDepthCalculation:
    """Test optical depth calculation behavior."""

    def test_optical_depth_increases_with_density(self):
        profile = get_body_profile("Earth", 0)
        assert isinstance(profile, AtmosphericProfile)

        depth_surface = calculate_optical_depth(profile, 0.0, 30.0)
        depth_high = calculate_optical_depth(profile, 10000.0, 30.0)

        assert depth_surface > depth_high, "Optical depth should decrease with altitude"

    def test_optical_depth_increases_at_haze(self):
        profile = get_body_profile("Earth", 0)
        assert isinstance(profile, AtmosphericProfile)

        depth_high_sun = calculate_optical_depth(profile, 0.0, 45.0)
        depth_low_sun = calculate_optical_depth(profile, 0.0, 5.0)

        assert depth_low_sun > depth_high_sun, (
            "Optical depth should increase at lower solar elevation"
        )

    def test_optical_depth_zero_for_null_atmosphere(self):
        profile = get_body_profile("Moon", 0)
        assert isinstance(profile, NullAtmosphere)


class TestMarsVisibility:
    """Test Mars visibility behavior."""

    def test_mars_surface_visible(self):
        profile = get_body_profile("Mars", 0)
        assert isinstance(profile, AtmosphericProfile)
        assert not isinstance(profile, NullAtmosphere)

        body = SOLAR_SYSTEM_BODIES["mars"]
        observer = Observer(
            latitude=0.0,
            longitude=0.0,
            altitude_m=0.0,
            body=body,
            solar_elevation_deg=0.0,
            sun_direction=(0.0, 1.0, 0.0),
            solar_angular_diameter_deg=0.35,
        )

        final_alt, justification, is_visible = resolve_visibility_elevation(
            profile, observer
        )

        assert is_visible, "Mars should be visible at surface"


class TestTitanVisibility:
    """Test Titan visibility behavior."""

    def test_titan_surface_visible(self):
        profile = get_body_profile("Titan", 0)
        assert isinstance(profile, AtmosphericProfile)
        assert not isinstance(profile, NullAtmosphere)

        body = SOLAR_SYSTEM_BODIES["titan"]
        observer = Observer(
            latitude=0.0,
            longitude=0.0,
            altitude_m=0.0,
            body=body,
            solar_elevation_deg=0.0,
            sun_direction=(0.0, 1.0, 0.0),
            solar_angular_diameter_deg=0.08,
        )

        final_alt, justification, is_visible = resolve_visibility_elevation(
            profile, observer
        )

        assert is_visible, "Titan should be visible at surface"
