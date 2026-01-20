import math
from dataclasses import dataclass

from ...models.atmosphere import AtmosphericProfile


def rayleigh_coefficient_co2(wavelength_nm: float) -> float:
    """Rayleigh scattering cross-section for CO₂ at standard conditions.
    Venus atmosphere is ~96.5% CO₂.

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
    """Combined Rayleigh coefficient for Venus' atmosphere.
    Weighted by surface mole fractions.
    """
    co2_fraction = 0.965
    n2_fraction = 0.035

    return co2_fraction * rayleigh_coefficient_co2(wavelength_nm)


def density_profile(altitude_m: float) -> float:
    """Venus' density profile using exponential barometric formula.
    Venus has an extremely dense atmosphere with smaller scale height.

    Args:
        altitude_m: Altitude in meters above reference surface

    Returns:
        Relative density (ρ/ρ₀) at given altitude
    """
    scale_height_m = 15800.0
    return math.exp(-altitude_m / scale_height_m)


def pressure_profile(altitude_m: float) -> float:
    """Venus' pressure profile using exponential barometric formula.

    Args:
        altitude_m: Altitude in meters above reference surface

    Returns:
        Relative pressure (P/P₀) at given altitude
    """
    scale_height_m = 15800.0
    return math.exp(-altitude_m / scale_height_m)


venus_profile = AtmosphericProfile(
    body_name="Venus",
    gas_composition={
        "CO2": 0.965,
        "N2": 0.035,
    },
    rayleigh_coefficient=combined_rayleigh_coefficient,
    mie_parameters={
        "cloud_drop_size_um": 2.0,
        "cloud_concentration": 2.0,
        "asymmetry_factor": 0.75,
    },
    density_profile=density_profile,
    pressure_profile=pressure_profile,
    refractive_index=1.000287,
    absorption_bands={
        "CO2": [
            (1950.0, 2100.0),
            (2700.0, 2950.0),
        ],
    },
)
