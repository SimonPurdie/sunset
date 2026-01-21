"""Tests for Optical Path Integrator."""

import math

import numpy as np
import pytest

from sunset.atmosphere.provider import get_body_profile
from sunset.models.atmosphere import NullAtmosphere
from sunset.optics.integrator import (
    compute_spectral_radiance,
    _compute_airmass,
    _compute_optical_depth,
    _rayleigh_phase_function,
    _henyey_greenstein_phase,
    _solar_irradiance,
)


class TestRayleighWavelengthScaling:
    """Test Rayleigh scattering scales approximately with 1/λ⁴."""

    def test_rayleigh_coefficient_scaling(self):
        """Rayleigh scattering coefficient follows 1/λ⁴ law."""
        profile = get_body_profile("Earth", 0)

        coeff_400 = profile.rayleigh_coefficient(400.0)
        coeff_550 = profile.rayleigh_coefficient(550.0)
        coeff_700 = profile.rayleigh_coefficient(700.0)

        ratio_400_550 = coeff_400 / coeff_550
        ratio_550_700 = coeff_550 / coeff_700

        expected_400_550 = (550.0 / 400.0) ** 4
        expected_550_700 = (700.0 / 550.0) ** 4

        assert math.isclose(ratio_400_550, expected_400_550, rel_tol=0.01), (
            f"Rayleigh coeff at 400nm vs 550nm: expected ~{expected_400_550:.2f}, got {ratio_400_550:.2f}"
        )
        assert math.isclose(ratio_550_700, expected_550_700, rel_tol=0.01), (
            f"Rayleigh coeff at 550nm vs 700nm: expected ~{expected_550_700:.2f}, got {ratio_550_700:.2f}"
        )


class TestShortWavelengthSkylight:
    """Test short wavelengths contribute more strongly to skylight at high solar elevations."""

    def test_blue_dominates_skylight_high_elevation(self):
        """Blue wavelengths dominate skylight at high solar elevation."""
        profile = get_body_profile("Earth", 0)

        wavelengths = np.arange(380, 781, 20.0)
        sun_direction = (0.0, 0.0, -1.0)
        body_radius_m = 6.371e6

        result = compute_spectral_radiance(
            profile, sun_direction, 0.0, 45.0, body_radius_m, wavelengths
        )

        skylight = result.skylight_contribution

        blue_wavelengths = wavelengths < 500
        red_wavelengths = wavelengths > 600

        avg_blue = np.mean(skylight[blue_wavelengths])
        avg_red = np.mean(skylight[red_wavelengths])

        assert avg_blue > avg_red, (
            f"Blue ({avg_blue:.2e}) should dominate skylight over red ({avg_red:.2e}) "
            "at high solar elevation"
        )


class TestOpticalDepthIncreasesHorizon:
    """Test optical depth increases as solar elevation approaches 0°."""

    def test_optical_depth_increases_near_horizon(self):
        """Optical depth increases as sun approaches horizon."""
        profile = get_body_profile("Earth", 0)

        wavelengths = np.array([550.0])
        sun_direction = (0.0, 0.0, -1.0)
        body_radius_m = 6.371e6

        result_high = compute_spectral_radiance(
            profile, sun_direction, 0.0, 45.0, body_radius_m, wavelengths
        )

        result_mid = compute_spectral_radiance(
            profile, sun_direction, 0.0, 10.0, body_radius_m, wavelengths
        )

        result_low = compute_spectral_radiance(
            profile, sun_direction, 0.0, 2.0, body_radius_m, wavelengths
        )

        direct_high = result_high.radiance[0]
        direct_mid = result_mid.radiance[0]
        direct_low = result_low.radiance[0]

        assert direct_low < direct_mid, (
            f"Direct radiance should decrease as sun approaches horizon: "
            f"low elevation {direct_low:.2e} < mid elevation {direct_mid:.2e}"
        )
        assert direct_mid < direct_high, (
            f"Direct radiance should decrease with lower elevation: "
            f"mid elevation {direct_mid:.2e} < high elevation {direct_high:.2e}"
        )

        airmass_high = _compute_airmass(math.radians(90.0 - 45.0), body_radius_m, 0.0)
        airmass_low = _compute_airmass(math.radians(90.0 - 2.0), body_radius_m, 0.0)

        assert airmass_low > airmass_high, (
            f"Airmass should increase near horizon: "
            f"low elevation {airmass_low:.2f} > high elevation {airmass_high:.2f}"
        )


class TestShortWavelengthAttenuation:
    """Test short wavelengths attenuate faster in thick atmospheres."""

    def test_blue_attenuates_more_in_thick_atmosphere(self):
        """Blue attenuates faster than red in thick atmospheres (e.g., Venus).

        Since blue scatters more strongly, the direct beam is attenuated more,
        so skylight proportion is higher for blue.
        """
        profile = get_body_profile("Venus", 0)

        wavelengths = np.array([400.0, 700.0])
        sun_direction = (0.0, 0.0, -1.0)
        body_radius_m = 6.052e6

        result = compute_spectral_radiance(
            profile, sun_direction, 0.0, 5.0, body_radius_m, wavelengths
        )

        blue_ratio = result.skylight_contribution[0] / result.radiance[0]
        red_ratio = result.skylight_contribution[1] / result.radiance[1]

        assert blue_ratio >= red_ratio, (
            f"Blue should have higher skylight ratio (more scattered) in thick atmosphere: "
            f"blue ratio {blue_ratio:.4f} >= red ratio {red_ratio:.4f}"
        )


class TestAirlessZeroSkylight:
    """Test airless bodies produce zero skylight."""

    def test_moon_zero_skylight(self):
        """Moon produces zero skylight."""
        profile = get_body_profile("Moon", 0)
        assert isinstance(profile, NullAtmosphere)

        wavelengths = np.arange(380, 781, 10.0)
        sun_direction = (0.0, 0.0, -1.0)
        body_radius_m = 1.737e6

        result = compute_spectral_radiance(
            profile, sun_direction, 0.0, 0.0, body_radius_m, wavelengths
        )

        assert np.all(result.skylight_contribution == 0), (
            "Skylight should be zero for airless bodies"
        )
        assert np.all(result.radiance == 0), (
            "Radiance should be zero for airless bodies (no atmosphere to scatter)"
        )

    def test_mercury_zero_skylight(self):
        """Mercury produces zero skylight."""
        profile = get_body_profile("Mercury", 0)
        assert isinstance(profile, NullAtmosphere)

        wavelengths = np.arange(380, 781, 10.0)
        sun_direction = (0.0, 0.0, -1.0)
        body_radius_m = 2.440e6

        result = compute_spectral_radiance(
            profile, sun_direction, 0.0, 0.0, body_radius_m, wavelengths
        )

        assert np.all(result.skylight_contribution == 0), (
            "Skylight should be zero for airless Mercury"
        )


class TestNoNegativeRadiance:
    """Test no negative radiance values."""

    def test_earth_no_negative_radiance(self):
        """Earth produces no negative radiance values."""
        profile = get_body_profile("Earth", 0)

        wavelengths = np.arange(380, 781, 10.0)
        sun_direction = (0.0, 0.0, -1.0)
        body_radius_m = 6.371e6

        result = compute_spectral_radiance(
            profile, sun_direction, 0.0, 30.0, body_radius_m, wavelengths
        )

        assert np.all(result.radiance >= 0), "Radiance should be non-negative"
        assert np.all(result.skylight_contribution >= 0), (
            "Skylight contribution should be non-negative"
        )

    def test_mars_no_negative_radiance(self):
        """Mars produces no negative radiance values."""
        profile = get_body_profile("Mars", 0)

        wavelengths = np.arange(380, 781, 10.0)
        sun_direction = (0.0, 0.0, -1.0)
        body_radius_m = 3.390e6

        result = compute_spectral_radiance(
            profile, sun_direction, 0.0, 10.0, body_radius_m, wavelengths
        )

        assert np.all(result.radiance >= 0), "Radiance should be non-negative"
        assert np.all(result.skylight_contribution >= 0), (
            "Skylight contribution should be non-negative"
        )

    def test_venus_no_negative_radiance(self):
        """Venus produces no negative radiance values."""
        profile = get_body_profile("Venus", 0)

        wavelengths = np.arange(380, 781, 10.0)
        sun_direction = (0.0, 0.0, -1.0)
        body_radius_m = 6.052e6

        result = compute_spectral_radiance(
            profile, sun_direction, 0.0, 5.0, body_radius_m, wavelengths
        )

        assert np.all(result.radiance >= 0), "Radiance should be non-negative"
        assert np.all(result.skylight_contribution >= 0), (
            "Skylight contribution should be non-negative"
        )


class TestRayleighPhaseFunction:
    """Test Rayleigh phase function behavior."""

    def test_rayleigh_phase_function_forward_backward(self):
        """Rayleigh phase is symmetric between forward and backward scattering."""
        phase_0 = _rayleigh_phase_function(0.0)
        phase_pi = _rayleigh_phase_function(math.pi)

        assert math.isclose(phase_0, phase_pi, rel_tol=1e-10), (
            "Rayleigh phase should be symmetric: P(0) = P(π)"
        )

    def test_rayleigh_phase_function_90_degrees(self):
        """Rayleigh phase at 90° should be exactly 3/(16π)."""
        phase_90 = _rayleigh_phase_function(math.pi / 2)
        expected = 3.0 / (16.0 * math.pi)

        assert math.isclose(phase_90, expected, rel_tol=1e-10), (
            f"Rayleigh phase at 90°: expected {expected:.6f}, got {phase_90:.6f}"
        )

    def test_rayleigh_phase_function_maximizes_at_0_180(self):
        """Rayleigh phase maximizes at 0° and 180°."""
        phase_0 = _rayleigh_phase_function(0.0)
        phase_90 = _rayleigh_phase_function(math.pi / 2)
        phase_45 = _rayleigh_phase_function(math.pi / 4)

        assert phase_0 > phase_90, "P(0°) > P(90°)"
        assert phase_0 > phase_45, "P(0°) > P(45°)"


class TestHenyeyGreensteinPhase:
    """Test Henyey-Greenstein phase function behavior."""

    def test_henyey_greenstein_symmetric_at_g_zero(self):
        """HG phase is symmetric when g=0."""
        phase_0 = _henyey_greenstein_phase(0.0, 0.0)
        phase_pi = _henyey_greenstein_phase(math.pi, 0.0)

        assert math.isclose(phase_0, phase_pi, rel_tol=1e-10), (
            "HG phase with g=0 should be symmetric"
        )

    def test_henyey_greenstein_forward_scattering(self):
        """HG phase with positive g peaks in forward direction."""
        phase_forward = _henyey_greenstein_phase(0.0, 0.7)
        phase_backward = _henyey_greenstein_phase(math.pi, 0.7)

        assert phase_forward > phase_backward, (
            f"Positive g should peak forward: {phase_forward:.6f} > {phase_backward:.6f}"
        )

    def test_henyey_greenstein_backward_scattering(self):
        """HG phase with negative g peaks in backward direction."""
        phase_forward = _henyey_greenstein_phase(0.0, -0.5)
        phase_backward = _henyey_greenstein_phase(math.pi, -0.5)

        assert phase_backward > phase_forward, (
            f"Negative g should peak backward: {phase_backward:.6f} > {phase_forward:.6f}"
        )


class TestSolarIrradiance:
    """Test solar irradiance model."""

    def test_solar_irradiance_positive(self):
        """Solar irradiance should be positive at all wavelengths."""
        wavelengths = np.arange(380, 781, 10.0)
        irradiance = _solar_irradiance(wavelengths)

        assert np.all(irradiance > 0), "Solar irradiance should be positive"

    def test_solar_irradiance_peaks_visible(self):
        """Solar irradiance should peak in visible range."""
        wavelengths = np.arange(200, 3000, 10.0)
        irradiance = _solar_irradiance(wavelengths)

        max_idx = np.argmax(irradiance)
        peak_wavelength = wavelengths[max_idx]

        assert 400 < peak_wavelength < 700, (
            f"Solar irradiance should peak in visible range, "
            f"but peak at {peak_wavelength:.0f}nm"
        )

    def test_solar_irradiance_decreases_at_uv_and_ir(self):
        """Solar irradiance should decrease away from visible."""
        wavelengths = np.array([300.0, 550.0, 2000.0])
        irradiance = _solar_irradiance(wavelengths)

        assert irradiance[0] < irradiance[1], "UV < visible peak"
        assert irradiance[2] < irradiance[1], "IR < visible peak"


class TestOpticalDepthIntegration:
    """Test optical depth computation."""

    def test_optical_depth_positive(self):
        """Optical depth should be positive."""
        profile = get_body_profile("Earth", 0)

        rayleigh_coeff = profile.rayleigh_coefficient(550.0)
        optical_depth = _compute_optical_depth(profile, rayleigh_coeff, 0.0, 1.0, 100)

        assert optical_depth > 0, "Optical depth should be positive"

    def test_optical_depth_increases_with_altitude_path(self):
        """Optical depth increases with longer atmospheric path."""
        profile = get_body_profile("Earth", 0)

        rayleigh_coeff = profile.rayleigh_coefficient(550.0)

        od_low = _compute_optical_depth(profile, rayleigh_coeff, 1000.0, 1.0, 100)
        od_high = _compute_optical_depth(profile, rayleigh_coeff, 1000.0, 3.0, 100)

        assert od_high > od_low, (
            f"Higher airmass should give higher optical depth: {od_high:.4f} > {od_low:.4f}"
        )


class TestAirmassComputation:
    """Test airmass computation."""

    def test_airmass_zenith(self):
        """Airmass at zenith should be approximately 1.0."""
        airmass = _compute_airmass(0.0, 6.371e6, 0.0)

        assert 0.9 < airmass < 1.2, (
            f"Airmass at zenith should be ~1.0, got {airmass:.2f}"
        )

    def test_airmass_increases_with_zenith_angle(self):
        """Airmass increases with zenith angle."""
        airmass_zenith = _compute_airmass(0.0, 6.371e6, 0.0)
        airmass_45 = _compute_airmass(math.radians(45.0), 6.371e6, 0.0)
        airmass_70 = _compute_airmass(math.radians(70.0), 6.371e6, 0.0)

        assert airmass_45 > airmass_zenith, (
            f"Airmass at 45° > zenith: {airmass_45:.2f} > {airmass_zenith:.2f}"
        )
        assert airmass_70 > airmass_45, (
            f"Airmass at 70° > 45°: {airmass_70:.2f} > {airmass_45:.2f}"
        )
