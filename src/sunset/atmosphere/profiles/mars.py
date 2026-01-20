import math
from dataclasses import dataclass

from ...models.atmosphere import AtmosphericProfile


def rayleigh_coefficient_co2(wavelength_nm: float) -> float:
    """Rayleigh scattering cross-section for CO₂ at standard conditions.
    Mars atmosphere is ~95% CO₂.

    Args:
        wavelength_nm: Wavelength in nanometers

    Returns:
        Scattering coefficient in m⁻¹ at standard conditions
    """
    wavelength_m = wavelength_nm * 1e-9
    ref_wavelength = 550e-9
    ref_coefficient = 8.7e-6

    return ref_coefficient * (ref_wavelength / wavelength_m) ** 4


def combined_rayleigh_coefficient(wavelength_nm: float) -> float:
    """Combined Rayleigh coefficient for Mars' atmosphere.
    Weighted by surface mole fractions.
    """
    co2_fraction = 0.953
    n2_fraction = 0.027
    ar_fraction = 0.016

    return co2_fraction * rayleigh_coefficient_co2(wavelength_nm)


def density_profile(altitude_m: float) -> float:
    """Mars' density profile using exponential barometric formula.
    Mars has a much thinner atmosphere with smaller scale height.

    Args:
        altitude_m: Altitude in meters above reference surface

    Returns:
        Relative density (ρ/ρ₀) at given altitude
    """
    scale_height_m = 11100.0
    return math.exp(-altitude_m / scale_height_m)


def pressure_profile(altitude_m: float) -> float:
    """Mars' pressure profile using exponential barometric formula.

    Args:
        altitude_m: Altitude in meters above reference surface

    Returns:
        Relative pressure (P/P₀) at given altitude
    """
    scale_height_m = 11100.0
    return math.exp(-altitude_m / scale_height_m)


mars_profile = AtmosphericProfile(
    body_name="Mars",
    gas_composition={
        "CO2": 0.957,
        "N2": 0.027,
        "Ar": 0.016,
    },
    rayleigh_coefficient=combined_rayleigh_coefficient,
    mie_parameters={
        "dust_size_um": 1.5,
        "dust_concentration": 0.5,
        "asymmetry_factor": 0.85,
    },
    density_profile=density_profile,
    pressure_profile=pressure_profile,
    refractive_index=1.000276,
    absorption_bands={},
)
