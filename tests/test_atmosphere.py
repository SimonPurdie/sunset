import math

import numpy as np
import pytest

from sunset.atmosphere.provider import get_body_profile
from sunset.models.atmosphere import AtmosphericProfile, NullAtmosphere


class TestDensityMonotonic:
    """Test that density decreases monotonically with altitude for all bodies."""

    def test_earth_density_monotonic(self):
        profile = get_body_profile("Earth", 0)
        assert isinstance(profile, AtmosphericProfile)

        altitudes = np.linspace(0, 50000, 100)
        densities = [profile.density_profile(alt) for alt in altitudes]

        for i in range(1, len(densities)):
            assert densities[i] <= densities[i - 1], (
                "Density should decrease with altitude"
            )

    def test_mars_density_monotonic(self):
        profile = get_body_profile("Mars", 0)
        assert isinstance(profile, AtmosphericProfile)

        altitudes = np.linspace(0, 50000, 100)
        densities = [profile.density_profile(alt) for alt in altitudes]

        for i in range(1, len(densities)):
            assert densities[i] <= densities[i - 1], (
                "Density should decrease with altitude"
            )

    def test_venus_density_monotonic(self):
        profile = get_body_profile("Venus", 0)
        assert isinstance(profile, AtmosphericProfile)

        altitudes = np.linspace(0, 50000, 100)
        densities = [profile.density_profile(alt) for alt in altitudes]

        for i in range(1, len(densities)):
            assert densities[i] <= densities[i - 1], (
                "Density should decrease with altitude"
            )

    def test_titan_density_monotonic(self):
        profile = get_body_profile("Titan", 0)
        assert isinstance(profile, AtmosphericProfile)

        altitudes = np.linspace(0, 50000, 100)
        densities = [profile.density_profile(alt) for alt in altitudes]

        for i in range(1, len(densities)):
            assert densities[i] <= densities[i - 1], (
                "Density should decrease with altitude"
            )


class TestPressureMonotonic:
    """Test that pressure decreases monotonically with altitude for all bodies."""

    def test_earth_pressure_monotonic(self):
        profile = get_body_profile("Earth", 0)
        assert isinstance(profile, AtmosphericProfile)

        altitudes = np.linspace(0, 50000, 100)
        pressures = [profile.pressure_profile(alt) for alt in altitudes]

        for i in range(1, len(pressures)):
            assert pressures[i] <= pressures[i - 1], (
                "Pressure should decrease with altitude"
            )

    def test_mars_pressure_monotonic(self):
        profile = get_body_profile("Mars", 0)
        assert isinstance(profile, AtmosphericProfile)

        altitudes = np.linspace(0, 50000, 100)
        pressures = [profile.pressure_profile(alt) for alt in altitudes]

        for i in range(1, len(pressures)):
            assert pressures[i] <= pressures[i - 1], (
                "Pressure should decrease with altitude"
            )

    def test_venus_pressure_monotonic(self):
        profile = get_body_profile("Venus", 0)
        assert isinstance(profile, AtmosphericProfile)

        altitudes = np.linspace(0, 50000, 100)
        pressures = [profile.pressure_profile(alt) for alt in altitudes]

        for i in range(1, len(pressures)):
            assert pressures[i] <= pressures[i - 1], (
                "Pressure should decrease with altitude"
            )

    def test_titan_pressure_monotonic(self):
        profile = get_body_profile("Titan", 0)
        assert isinstance(profile, AtmosphericProfile)

        altitudes = np.linspace(0, 50000, 100)
        pressures = [profile.pressure_profile(alt) for alt in altitudes]

        for i in range(1, len(pressures)):
            assert pressures[i] <= pressures[i - 1], (
                "Pressure should decrease with altitude"
            )


class TestEarthMatchesPublished:
    """Test Earth profile matches published atmospheric ranges."""

    def test_gas_composition_sum_to_one(self):
        profile = get_body_profile("Earth", 0)
        assert isinstance(profile, AtmosphericProfile)

        total = sum(profile.gas_composition.values())
        assert math.isclose(total, 1.0, rel_tol=1e-10), (
            "Gas composition should sum to 1.0"
        )

    def test_primary_gases_present(self):
        profile = get_body_profile("Earth", 0)
        assert isinstance(profile, AtmosphericProfile)

        assert "N2" in profile.gas_composition
        assert "O2" in profile.gas_composition
        assert "Ar" in profile.gas_composition
        assert "CO2" in profile.gas_composition

    def test_gas_composition_within_published_ranges(self):
        profile = get_body_profile("Earth", 0)
        assert isinstance(profile, AtmosphericProfile)

        n2_range = (0.780, 0.784)
        o2_range = (0.209, 0.210)
        ar_range = (0.009, 0.0095)

        assert n2_range[0] <= profile.gas_composition["N2"] <= n2_range[1]
        assert o2_range[0] <= profile.gas_composition["O2"] <= o2_range[1]
        assert ar_range[0] <= profile.gas_composition["Ar"] <= ar_range[1]

    def test_scale_height_reasonable(self):
        profile = get_body_profile("Earth", 0)
        assert isinstance(profile, AtmosphericProfile)

        density_surface = profile.density_profile(0)
        density_10km = profile.density_profile(10000)
        density_20km = profile.density_profile(20000)

        ratio_10km = density_10km / density_surface
        ratio_20km = density_20km / density_surface

        expected_10km = math.exp(-10000 / 8500)
        expected_20km = math.exp(-20000 / 8500)

        assert math.isclose(ratio_10km, expected_10km, rel_tol=0.01)
        assert math.isclose(ratio_20km, expected_20km, rel_tol=0.01)


class TestMarsMatchesPublished:
    """Test Mars profile matches published atmospheric ranges."""

    def test_gas_composition_sum_to_one(self):
        profile = get_body_profile("Mars", 0)
        assert isinstance(profile, AtmosphericProfile)

        total = sum(profile.gas_composition.values())
        assert math.isclose(total, 1.0, rel_tol=1e-10)

    def test_primary_gas_is_co2(self):
        profile = get_body_profile("Mars", 0)
        assert isinstance(profile, AtmosphericProfile)

        assert profile.gas_composition["CO2"] > 0.95
        assert "N2" in profile.gas_composition
        assert "Ar" in profile.gas_composition

    def test_scale_height_reasonable(self):
        profile = get_body_profile("Mars", 0)
        assert isinstance(profile, AtmosphericProfile)

        density_surface = profile.density_profile(0)
        density_10km = profile.density_profile(10000)

        ratio_10km = density_10km / density_surface
        expected_10km = math.exp(-10000 / 11100)

        assert math.isclose(ratio_10km, expected_10km, rel_tol=0.01)


class TestVenusMatchesPublished:
    """Test Venus profile matches published atmospheric ranges."""

    def test_gas_composition_sum_to_one(self):
        profile = get_body_profile("Venus", 0)
        assert isinstance(profile, AtmosphericProfile)

        total = sum(profile.gas_composition.values())
        assert math.isclose(total, 1.0, rel_tol=1e-10)

    def test_primary_gas_is_co2(self):
        profile = get_body_profile("Venus", 0)
        assert isinstance(profile, AtmosphericProfile)

        assert profile.gas_composition["CO2"] > 0.96
        assert "N2" in profile.gas_composition

    def test_scale_height_reasonable(self):
        profile = get_body_profile("Venus", 0)
        assert isinstance(profile, AtmosphericProfile)

        density_surface = profile.density_profile(0)
        density_10km = profile.density_profile(10000)

        ratio_10km = density_10km / density_surface
        expected_10km = math.exp(-10000 / 15800)

        assert math.isclose(ratio_10km, expected_10km, rel_tol=0.01)


class TestTitanMatchesPublished:
    """Test Titan profile matches published atmospheric ranges."""

    def test_gas_composition_sum_to_one(self):
        profile = get_body_profile("Titan", 0)
        assert isinstance(profile, AtmosphericProfile)

        total = sum(profile.gas_composition.values())
        assert math.isclose(total, 1.0, rel_tol=1e-10)

    def test_primary_gas_is_n2(self):
        profile = get_body_profile("Titan", 0)
        assert isinstance(profile, AtmosphericProfile)

        assert profile.gas_composition["N2"] > 0.98
        assert "CH4" in profile.gas_composition

    def test_scale_height_reasonable(self):
        profile = get_body_profile("Titan", 0)
        assert isinstance(profile, AtmosphericProfile)

        density_surface = profile.density_profile(0)
        density_10km = profile.density_profile(10000)

        ratio_10km = density_10km / density_surface
        expected_10km = math.exp(-10000 / 21000)

        assert math.isclose(ratio_10km, expected_10km, rel_tol=0.01)


class TestAirlessNull:
    """Test that airless bodies return NullAtmosphere."""

    def test_mercury_returns_null_atmosphere(self):
        profile = get_body_profile("Mercury", 0)
        assert isinstance(profile, NullAtmosphere)
        assert profile.body_name == "Mercury"

    def test_moon_returns_null_atmosphere(self):
        profile = get_body_profile("Moon", 0)
        assert isinstance(profile, NullAtmosphere)
        assert profile.body_name == "Moon"


class TestRayleighScattering:
    """Test Rayleigh scattering coefficients follow λ⁻⁴ law."""

    def test_earth_rayleigh_scaling(self):
        profile = get_body_profile("Earth", 0)
        assert isinstance(profile, AtmosphericProfile)

        wavelength_400 = profile.rayleigh_coefficient(400)
        wavelength_800 = profile.rayleigh_coefficient(800)

        ratio = wavelength_400 / wavelength_800
        expected_ratio = 16.0  # (800/400)⁴ = 2⁴ = 16

        assert math.isclose(ratio, expected_ratio, rel_tol=0.1)

    def test_mars_rayleigh_scaling(self):
        profile = get_body_profile("Mars", 0)
        assert isinstance(profile, AtmosphericProfile)

        wavelength_400 = profile.rayleigh_coefficient(400)
        wavelength_800 = profile.rayleigh_coefficient(800)

        ratio = wavelength_400 / wavelength_800
        expected_ratio = 16.0

        assert math.isclose(ratio, expected_ratio, rel_tol=0.1)

    def test_rayleigh_positive(self):
        profile = get_body_profile("Earth", 0)
        assert isinstance(profile, AtmosphericProfile)

        for wavelength in [400, 500, 600, 700]:
            coeff = profile.rayleigh_coefficient(wavelength)
            assert coeff > 0


class TestRefractiveIndex:
    """Test refractive indices are physically reasonable."""

    def test_earth_refractive_index(self):
        profile = get_body_profile("Earth", 0)
        assert isinstance(profile, AtmosphericProfile)

        assert 1.0 < profile.refractive_index < 1.01

    def test_mars_refractive_index(self):
        profile = get_body_profile("Mars", 0)
        assert isinstance(profile, AtmosphericProfile)

        assert 1.0 < profile.refractive_index < 1.01

    def test_venus_refractive_index(self):
        profile = get_body_profile("Venus", 0)
        assert isinstance(profile, AtmosphericProfile)

        assert 1.0 < profile.refractive_index < 1.01

    def test_titan_refractive_index(self):
        profile = get_body_profile("Titan", 0)
        assert isinstance(profile, AtmosphericProfile)

        assert 1.0 < profile.refractive_index < 1.01
