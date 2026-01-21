"""Exposure and tone mapping for converting linear colors to displayable range."""

import numpy as np


def apply_exposure_and_tonemap(xyz_linear: np.ndarray) -> np.ndarray:
    """Apply exposure and tone mapping to linear XYZ values.

    Uses automatic exposure based on luminance, then applies
    Reinhard tone mapping to preserve hue.

    Args:
        xyz_linear: Linear XYZ tristimulus values

    Returns:
        np.ndarray: Tone-mapped XYZ values in [0, 1] range
    """
    Y = xyz_linear[1]

    if Y <= 0:
        return np.zeros(3)

    exposure_factor = _compute_auto_exposure(Y)

    xyz_exposed = xyz_linear * exposure_factor

    xyz_tonemapped = _reinhard_tonemap(xyz_exposed)

    return xyz_tonemapped


def _compute_auto_exposure(Y: float) -> float:
    """Compute automatic exposure factor based on luminance.

    Targets middle gray (Y = 0.18) for the scene.

    Args:
        Y: Luminance value

    Returns:
        float: Exposure factor
    """
    target_middle_gray = 0.18
    if Y <= 0:
        return 1.0

    return target_middle_gray / Y


def _reinhard_tonemap(xyz: np.ndarray) -> np.ndarray:
    """Apply Reinhard tone mapping operator.

    Reinhard: f(x) = x / (1 + x)
    Preserves hue by applying same mapping to all channels.

    Args:
        xyz: Linear XYZ values (may be > 1.0)

    Returns:
        np.ndarray: Tone-mapped values in [0, 1] range
    """
    return xyz / (1.0 + xyz)


def apply_simple_exposure(
    xyz_linear: np.ndarray, exposure_stops: float = 0.0
) -> np.ndarray:
    """Apply simple exposure adjustment without tone mapping.

    Args:
        xyz_linear: Linear XYZ tristimulus values
        exposure_stops: Exposure adjustment in stops (log2 scale)

    Returns:
        np.ndarray: Exposed linear XYZ values
    """
    exposure_factor = 2.0**exposure_stops
    return xyz_linear * exposure_factor
