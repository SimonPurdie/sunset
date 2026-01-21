"""XYZ to sRGB color encoding with gamma correction."""

import numpy as np


# XYZ to linear sRGB transformation matrix (D65 illuminant)
_XYZ_TO_SRGB_MATRIX = np.array(
    [
        [3.2406, -1.5372, -0.4986],
        [-0.9689, 1.8758, 0.0415],
        [0.0557, -0.2040, 1.0570],
    ]
)


def xyz_to_srgb(xyz: np.ndarray) -> np.ndarray:
    """Convert XYZ to sRGB with gamma encoding.

    Args:
        xyz: XYZ tristimulus values (should be in [0, 1] range after tone mapping)

    Returns:
        np.ndarray: sRGB values in [0, 1] range
    """
    linear_srgb = _xyz_to_linear_srgb(xyz)

    srgb_encoded = _apply_srgb_gamma(linear_srgb)

    return srgb_encoded


def _xyz_to_linear_srgb(xyz: np.ndarray) -> np.ndarray:
    """Convert XYZ to linear sRGB using transformation matrix.

    Args:
        xyz: XYZ tristimulus values

    Returns:
        np.ndarray: Linear sRGB values
    """
    return _XYZ_TO_SRGB_MATRIX @ xyz


def _apply_srgb_gamma(linear_srgb: np.ndarray) -> np.ndarray:
    """Apply sRGB gamma encoding.

    sRGB transfer function (IEC 61966-2-1:1999):
    - If linear_sRGB <= 0.0031308: 12.92 * linear_sRGB
    - Otherwise: 1.055 * linear_sRGB^(1/2.4) - 0.055

    Args:
        linear_srgb: Linear sRGB values in [0, 1] range

    Returns:
        np.ndarray: Gamma-encoded sRGB values
    """
    threshold = 0.0031308
    a = 12.92
    b = 1.055
    c = 1.0 / 2.4
    d = 0.055

    result = np.where(
        linear_srgb <= threshold, a * linear_srgb, b * np.power(linear_srgb, c) - d
    )

    return np.clip(result, 0.0, 1.0)


def linear_to_srgb_16bit(srgb_normalized: np.ndarray) -> np.ndarray:
    """Convert normalized sRGB values to 16-bit integer range.

    Args:
        srgb_normalized: sRGB values in [0, 1] range

    Returns:
        np.ndarray: sRGB values scaled to 16-bit range [0, 65535]
    """
    return np.clip(srgb_normalized * 65535.0 + 0.5, 0, 65535).astype(np.uint16)


def linear_to_srgb_8bit(srgb_normalized: np.ndarray) -> np.ndarray:
    """Convert normalized sRGB values to 8-bit integer range.

    Args:
        srgb_normalized: sRGB values in [0, 1] range

    Returns:
        np.ndarray: sRGB values scaled to 8-bit range [0, 255]
    """
    return np.clip(srgb_normalized * 255.0 + 0.5, 0, 255).astype(np.uint8)
