import math
from dataclasses import dataclass

from ...models.atmosphere import AtmosphericProfile


def rayleigh_coefficient_n2(wavelength_nm: float) -> float:
    """Rayleigh scattering cross-section for N₂ at standard conditions.
    Titan atmosphere is ~98.4% N₂.

    Args:
        wavelength_nm: Wavelength in nanometers

    Returns:
        Scattering coefficient in m⁻¹ at standard conditions
    """
    wavelength_m = wavelength_nm * 1e-9
    ref_wavelength = 550e-9
    ref_coefficient = 1.54e-5

    return ref_coefficient * (ref_wavelength / wavelength_m) ** 4


def rayleigh_coefficient_ch4(wavelength_nm: float) -> float:
    """Rayleigh scattering cross-section for CH₄ at standard conditions.

    Args:
        wavelength_nm: Wavelength in nanometers

    Returns:
        Scattering coefficient in m⁻¹ at standard conditions
    """
    wavelength_m = wavelength_nm * 1e-9
    ref_wavelength = 550e-9
    ref_coefficient = 7.4e-6

    return ref_coefficient * (ref_wavelength / wavelength_m) ** 4


def combined_rayleigh_coefficient(wavelength_nm: float) -> float:
    """Combined Rayleigh coefficient for Titan's atmosphere.
    Weighted by surface mole fractions.
    """
    n2_fraction = 0.984
    ch4_fraction = 0.014

    return n2_fraction * rayleigh_coefficient_n2(
        wavelength_nm
    ) + ch4_fraction * rayleigh_coefficient_ch4(wavelength_nm)


def density_profile(altitude_m: float) -> float:
    """Titan's density profile using exponential barometric formula.
    Titan has a thick atmosphere with large scale height due to low gravity.

    Args:
        altitude_m: Altitude in meters above reference surface

    Returns:
        Relative density (ρ/ρ₀) at given altitude
    """
    scale_height_m = 21000.0
    return math.exp(-altitude_m / scale_height_m)


def pressure_profile(altitude_m: float) -> float:
    """Titan's pressure profile using exponential barometric formula.

    Args:
        altitude_m: Altitude in meters above reference surface

    Returns:
        Relative pressure (P/P₀) at given altitude
    """
    scale_height_m = 21000.0
    return math.exp(-altitude_m / scale_height_m)


titan_profile = AtmosphericProfile(
    body_name="Titan",
    gas_composition={
        "N2": 0.984,
        "CH4": 0.014,
        "Ar": 0.002,
    },
    rayleigh_coefficient=combined_rayleigh_coefficient,
    mie_parameters={
        "haze_size_um": 0.1,
        "haze_concentration": 3.0,
        "asymmetry_factor": 0.8,
    },
    density_profile=density_profile,
    pressure_profile=pressure_profile,
    refractive_index=1.000312,
    absorption_bands={
        "CH4": [
            (720.0, 950.0),
            (1100.0, 1350.0),
            (1650.0, 1750.0),
        ],
    },
)
