import math
from dataclasses import dataclass

from ...models.atmosphere import AtmosphericProfile


def rayleigh_coefficient_n2(wavelength_nm: float) -> float:
    """Rayleigh scattering cross-section for N₂ at standard conditions.
    Follows λ⁻⁴ wavelength dependence.

    Args:
        wavelength_nm: Wavelength in nanometers

    Returns:
        Scattering coefficient in m⁻¹ at standard conditions
    """
    wavelength_m = wavelength_nm * 1e-9
    ref_wavelength = 550e-9
    ref_coefficient = 1.54e-5

    return ref_coefficient * (ref_wavelength / wavelength_m) ** 4


def rayleigh_coefficient_o2(wavelength_nm: float) -> float:
    """Rayleigh scattering cross-section for O₂ at standard conditions.
    Follows λ⁻⁴ wavelength dependence.
    """
    wavelength_m = wavelength_nm * 1e-9
    ref_wavelength = 550e-9
    ref_coefficient = 1.19e-5

    return ref_coefficient * (ref_wavelength / wavelength_m) ** 4


def rayleigh_coefficient_ar(wavelength_nm: float) -> float:
    """Rayleigh scattering cross-section for Ar at standard conditions.
    Follows λ⁻⁴ wavelength dependence.
    """
    wavelength_m = wavelength_nm * 1e-9
    ref_wavelength = 550e-9
    ref_coefficient = 7.1e-6

    return ref_coefficient * (ref_wavelength / wavelength_m) ** 4


def rayleigh_coefficient_co2(wavelength_nm: float) -> float:
    """Rayleigh scattering cross-section for CO₂ at standard conditions.
    Follows λ⁻⁴ wavelength dependence.
    """
    wavelength_m = wavelength_nm * 1e-9
    ref_wavelength = 550e-9
    ref_coefficient = 8.7e-6

    return ref_coefficient * (ref_wavelength / wavelength_m) ** 4


def combined_rayleigh_coefficient(wavelength_nm: float) -> float:
    """Combined Rayleigh coefficient for Earth's atmosphere (N₂ + O₂ + Ar + CO₂).
    Weighted by sea-level mole fractions.
    """
    n2_fraction = 0.7808
    o2_fraction = 0.2095
    ar_fraction = 0.0093
    co2_fraction = 0.0004

    return (
        n2_fraction * rayleigh_coefficient_n2(wavelength_nm)
        + o2_fraction * rayleigh_coefficient_o2(wavelength_nm)
        + ar_fraction * rayleigh_coefficient_ar(wavelength_nm)
        + co2_fraction * rayleigh_coefficient_co2(wavelength_nm)
    )


def density_profile(altitude_m: float) -> float:
    """Earth's density profile using exponential barometric formula.

    Args:
        altitude_m: Altitude in meters above sea level

    Returns:
        Relative density (ρ/ρ₀) at given altitude
    """
    scale_height_m = 8500.0
    return math.exp(-altitude_m / scale_height_m)


def pressure_profile(altitude_m: float) -> float:
    """Earth's pressure profile using exponential barometric formula.

    Args:
        altitude_m: Altitude in meters above sea level

    Returns:
        Relative pressure (P/P₀) at given altitude
    """
    scale_height_m = 8500.0
    return math.exp(-altitude_m / scale_height_m)


earth_profile = AtmosphericProfile(
    body_name="Earth",
    gas_composition={
        "N2": 0.7808,
        "O2": 0.2095,
        "Ar": 0.0093,
        "CO2": 0.0004,
    },
    rayleigh_coefficient=combined_rayleigh_coefficient,
    mie_parameters={
        "aerosol_size_um": 0.5,
        "aerosol_concentration": 1.0,
        "asymmetry_factor": 0.7,
    },
    density_profile=density_profile,
    pressure_profile=pressure_profile,
    refractive_index=1.000293,
    absorption_bands={
        "O2": [
            (760.0, 770.0),
            (687.0, 692.0),
        ],
        "H2O": [
            (935.0, 985.0),
            (1125.0, 1165.0),
        ],
    },
)
