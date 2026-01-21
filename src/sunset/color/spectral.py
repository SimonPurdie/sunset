"""Spectral to XYZ color conversion using CIE standard observer."""

import numpy as np
from colour import SpectralDistribution, sd_to_XYZ

from ..models.spectral import SpectralRadiance


def spectral_to_xyz(spectral_radiance: SpectralRadiance) -> np.ndarray:
    """Convert spectral radiance to CIE XYZ tristimulus values.

    Uses CIE 1931 2° standard observer color matching functions.
    Integrates spectral radiance × CMF over wavelength.

    Args:
        spectral_radiance: Spectral radiance data with wavelengths and radiance

    Returns:
        np.ndarray: XYZ tristimulus values [X, Y, Z]
    """
    wavelengths = spectral_radiance.wavelengths
    radiance = spectral_radiance.radiance

    sd = SpectralDistribution(
        data=radiance,
        domain=wavelengths,
        name="spectral_radiance",
    )

    xyz_result = sd_to_XYZ(sd, method="ASTM E308")

    return np.array([xyz_result[0], xyz_result[1], xyz_result[2]])
