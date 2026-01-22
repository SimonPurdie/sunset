"""Star catalog module providing bright star positions and magnitudes."""

import numpy as np

# Bright stars catalog (RA in hours, Dec in degrees, V magnitude)
# RA is right ascension (0-24 hours), Dec is declination (-90 to +90 degrees)
# Selected brightest and most recognizable stars
BRIGHT_STARS = np.array(
    [
        (6.75, -16.72, -1.46),  # Sirius
        (14.26, -60.83, -0.74),  # Canopus
        (5.92, 7.41, 0.03),  # Arcturus
        (13.97, -11.15, -0.05),  # Vega
        (5.55, 7.24, 0.12),  # Capella
        (6.75, -16.72, -1.46),  # Rigel
        (5.24, 6.35, 0.18),  # Procyon
        (5.92, 7.41, 0.03),  # Betelgeuse
        (13.70, 49.31, 0.40),  # Altair
        (18.62, 38.78, 0.77),  # Deneb
        (20.69, 45.28, 0.85),  # Alderamin
        (5.92, 7.41, 0.03),  # Aldebaran
        (13.25, -43.20, 0.60),  # Acrux
        (5.92, 7.41, 0.03),  # Antares
        (11.06, 61.75, 0.08),  # Pollux
        (10.08, -11.97, 0.98),  # Fomalhaut
        (16.71, -26.30, 0.42),  # Deneb
        (14.15, 19.18, 0.35),  # Spica
        (7.75, 28.03, 0.50),  # Regulus
        (20.41, 45.28, 1.25),  # Enif
    ],
    dtype=np.float32,
)


def get_bright_stars(max_magnitude: float = 2.0) -> np.ndarray:
    """Get bright stars from the catalog.

    Args:
        max_magnitude: Maximum visual magnitude (lower = brighter)

    Returns:
        Array of (RA_hours, Dec_deg, magnitude) for stars brighter than max_magnitude
    """
    mask = BRIGHT_STARS[:, 2] <= max_magnitude
    return BRIGHT_STARS[mask]


def ra_hours_to_radians(ra_hours: np.ndarray) -> np.ndarray:
    """Convert right ascension from hours to radians.

    Args:
        ra_hours: Right ascension in hours (0-24)

    Returns:
        Right ascension in radians (0-2π)
    """
    return ra_hours * (2.0 * np.pi / 24.0)


def dec_deg_to_radians(dec_deg: np.ndarray) -> np.ndarray:
    """Convert declination from degrees to radians.

    Args:
        dec_deg: Declination in degrees (-90 to +90)

    Returns:
        Declination in radians (-π/2 to π/2)
    """
    return np.radians(dec_deg)


def ra_dec_to_cartesian(ra_hours: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    """Convert RA/Dec to unit vector in Cartesian coordinates (ICRS frame).

    Args:
        ra_hours: Right ascension in hours
        dec_deg: Declination in degrees

    Returns:
        Array of unit vectors with shape (N, 3) in ICRS frame
    """
    ra = ra_hours_to_radians(ra_hours)
    dec = dec_deg_to_radians(dec_deg)

    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)

    return np.stack([x, y, z], axis=1)


def apparent_magnitude_to_brightness(magnitude: np.ndarray) -> np.ndarray:
    """Convert apparent magnitude to relative brightness.

    Brightness follows: B = 10^(-0.4 * magnitude)
    Lower magnitude = brighter star (e.g., magnitude -1.5 is very bright)

    Args:
        magnitude: Apparent V magnitude

    Returns:
        Relative brightness (higher = brighter)
    """
    return 10.0 ** (-0.4 * magnitude)
