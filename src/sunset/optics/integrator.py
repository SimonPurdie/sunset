"""Optical Path Integrator for computing spectral radiance."""

import math
from typing import Tuple

import numpy as np

from ..models.atmosphere import AtmosphericProfile, NullAtmosphere
from ..models.spectral import SpectralRadiance


def compute_spectral_radiance(
    atmospheric_profile: AtmosphericProfile | NullAtmosphere,
    sun_direction: Tuple[float, float, float],
    observer_altitude_m: float,
    solar_elevation_deg: float,
    body_radius_m: float,
    wavelengths_nm: np.ndarray | None = None,
    num_integration_steps: int = 100,
) -> SpectralRadiance:
    """Compute spectral radiance reaching observer.

    Args:
        atmospheric_profile: Atmospheric profile for the body
        sun_direction: Unit vector pointing toward sun (body-fixed coords)
        observer_altitude_m: Observer altitude in meters
        solar_elevation_deg: Solar elevation angle in degrees
        body_radius_m: Body radius in meters
        wavelengths_nm: Wavelength array in nm (default: 380-780nm in 10nm steps)
        num_integration_steps: Number of steps for path integration

    Returns:
        SpectralRadiance with wavelength-dependent radiance and skylight contribution
    """
    if wavelengths_nm is None:
        wavelengths_nm = np.arange(380, 781, 10.0)

    if isinstance(atmospheric_profile, NullAtmosphere):
        return _null_atmosphere_radiance(wavelengths_nm)

    return _atmospheric_radiance(
        atmospheric_profile,
        sun_direction,
        observer_altitude_m,
        solar_elevation_deg,
        body_radius_m,
        wavelengths_nm,
        num_integration_steps,
    )


def _null_atmosphere_radiance(wavelengths_nm: np.ndarray) -> SpectralRadiance:
    """Return zero radiance for airless bodies."""
    num_wavelengths = len(wavelengths_nm)
    return SpectralRadiance(
        wavelengths=wavelengths_nm,
        radiance=np.zeros(num_wavelengths),
        skylight_contribution=np.zeros(num_wavelengths),
    )


def _atmospheric_radiance(
    profile: AtmosphericProfile,
    sun_direction: Tuple[float, float, float],
    observer_altitude_m: float,
    solar_elevation_deg: float,
    body_radius_m: float,
    wavelengths_nm: np.ndarray,
    num_integration_steps: int,
) -> SpectralRadiance:
    """Compute radiance with atmospheric scattering."""
    num_wavelengths = len(wavelengths_nm)

    solar_zenith_angle_rad = math.radians(90.0 - solar_elevation_deg)

    airmass = _compute_airmass(
        solar_zenith_angle_rad, body_radius_m, observer_altitude_m
    )

    optical_depth = np.zeros(num_wavelengths)
    transmission = np.ones(num_wavelengths)
    scattering_optical_depth = np.zeros(num_wavelengths)

    for i, wavelength_nm in enumerate(wavelengths_nm):
        rayleigh_coeff = profile.rayleigh_coefficient(wavelength_nm)

        total_optical_depth = _compute_optical_depth(
            profile, rayleigh_coeff, observer_altitude_m, airmass, num_integration_steps
        )

        optical_depth[i] = total_optical_depth
        transmission[i] = math.exp(-total_optical_depth)
        scattering_optical_depth[i] = total_optical_depth * _scattering_fraction(
            profile
        )

    direct_radiance = transmission * _solar_irradiance(wavelengths_nm)

    skylight_contribution = _compute_skylight(
        wavelengths_nm,
        scattering_optical_depth,
        solar_zenith_angle_rad,
        profile,
    )

    total_radiance = direct_radiance + skylight_contribution

    return SpectralRadiance(
        wavelengths=wavelengths_nm,
        radiance=total_radiance,
        skylight_contribution=skylight_contribution,
    )


def _compute_airmass(
    solar_zenith_angle_rad: float,
    body_radius_m: float,
    observer_altitude_m: float,
) -> float:
    """Compute airmass factor for spherical atmosphere.

    Uses Kasten & Young formula modified for spherical geometry.
    """
    cos_zenith = math.cos(solar_zenith_angle_rad)

    altitude_fraction = observer_altitude_m / body_radius_m

    if cos_zenith > 0:
        airmass = 1.0 / (
            cos_zenith
            + 0.50572 * (96.07995 - solar_zenith_angle_rad * 180 / math.pi) ** -1.6364
        )
    else:
        sin_horizon = body_radius_m / (body_radius_m + observer_altitude_m)
        horizon_zenith = math.asin(sin_horizon)
        if solar_zenith_angle_rad > horizon_zenith:
            airmass = 40.0
        else:
            airmass = 1.0 / math.cos(solar_zenith_angle_rad)

    airmass = min(airmass, 40.0)
    return airmass * (1.0 - altitude_fraction * 0.5)


def _compute_optical_depth(
    profile: AtmosphericProfile,
    rayleigh_coeff: float,
    observer_altitude_m: float,
    airmass: float,
    num_steps: int,
) -> float:
    """Compute optical depth along sun path through atmosphere.

    Integrates scattering coefficient along line-of-sight to sun.
    """
    if profile.density_profile is None:
        density_fn = lambda alt: 1.0
    else:
        density_fn = profile.density_profile

    max_altitude_m = min(observer_altitude_m + 100000, 100000)

    altitudes = np.linspace(observer_altitude_m, max_altitude_m, num_steps)
    step_size = (max_altitude_m - observer_altitude_m) / (num_steps - 1)

    optical_depth = 0.0

    for altitude in altitudes:
        relative_density = density_fn(altitude)
        scattering_coeff = rayleigh_coeff * relative_density

        if profile.mie_parameters:
            aerosol_concentration = profile.mie_parameters.get(
                "aerosol_concentration", 1.0
            )
            aerosol_contribution = 5e-6 * aerosol_concentration * relative_density
            scattering_coeff += aerosol_contribution

        optical_depth += scattering_coeff * step_size

    return optical_depth * airmass


def _scattering_fraction(profile: AtmosphericProfile) -> float:
    """Return fraction of optical depth due to scattering vs absorption."""
    if not profile.absorption_bands:
        return 1.0
    return 0.95


def _solar_irradiance(wavelengths_nm: np.ndarray) -> np.ndarray:
    """Return solar spectral irradiance at top of atmosphere.

    Simplified model using solar constant with blackbody spectrum.
    """
    solar_constant = 1361.0

    temperature = 5778.0
    h = 6.62607015e-34
    c = 2.99792458e8
    k = 1.380649e-23

    wavelengths_m = wavelengths_nm * 1e-9

    spectral_radiance = (2 * h * c**2 / wavelengths_m**5) / (
        np.exp(h * c / (wavelengths_m * k * temperature)) - 1
    )

    if len(wavelengths_nm) == 1:
        return spectral_radiance * 1e-9

    integrated_irradiance = np.trapezoid(spectral_radiance, wavelengths_m)

    scale_factor = solar_constant / integrated_irradiance

    return spectral_radiance * scale_factor


def _compute_skylight(
    wavelengths_nm: np.ndarray,
    scattering_optical_depth: np.ndarray,
    solar_zenith_angle_rad: float,
    profile: AtmosphericProfile,
) -> np.ndarray:
    """Compute skylight contribution using single-scattering approximation."""
    single_scattering_albedo = _single_scattering_albedo(profile)

    phase_function = _rayleigh_phase_function(solar_zenith_angle_rad)

    skylight_factor = (
        single_scattering_albedo
        * phase_function
        * (1 - np.exp(-scattering_optical_depth))
    )

    skylight_contribution = skylight_factor * _solar_irradiance(wavelengths_nm) * 0.3

    if profile.mie_parameters:
        asymmetry_factor = profile.mie_parameters.get("asymmetry_factor", 0.7)
        mie_phase = _henyey_greenstein_phase(solar_zenith_angle_rad, asymmetry_factor)
        mie_skylight = 0.1 * mie_phase * _solar_irradiance(wavelengths_nm)
        skylight_contribution += mie_skylight

    return skylight_contribution


def _single_scattering_albedo(profile: AtmosphericProfile) -> float:
    """Return single-scattering albedo (scattering / total extinction)."""
    if not profile.absorption_bands:
        return 1.0
    return 0.98


def _rayleigh_phase_function(scattering_angle_rad: float) -> float:
    """Rayleigh phase function: P(θ) = (3/16π) * (1 + cos²θ)."""
    cos_theta = math.cos(scattering_angle_rad)
    return 3.0 / (16.0 * math.pi) * (1.0 + cos_theta**2)


def _henyey_greenstein_phase(
    scattering_angle_rad: float, asymmetry_factor: float
) -> float:
    """Henyey-Greenstein phase function for Mie scattering."""
    if abs(asymmetry_factor) < 1e-6:
        return 1.0 / (4.0 * math.pi)

    g = asymmetry_factor
    cos_theta = math.cos(scattering_angle_rad)

    denominator = (1 + g**2 - 2 * g * cos_theta) ** 1.5
    return (1 - g**2) / (4 * math.pi * denominator)
