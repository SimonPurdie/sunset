"""PNG metadata embedding and extraction for sunset scenes."""

from typing import Dict, Optional

import numpy as np
from PIL import Image, PngImagePlugin

from ..models.scene import Scene


def embed_metadata(image_array: np.ndarray, scene: Scene, output_path: str) -> None:
    """Save a 16-bit PNG image with embedded metadata.

    Args:
        image_array: 16-bit RGB image array with shape (height, width, 3)
        scene: Scene metadata object
        output_path: Path where to save the PNG file
    """
    if image_array.dtype != np.uint16:
        raise ValueError("Image array must be uint16 for 16-bit PNG output")

    if len(image_array.shape) != 3 or image_array.shape[2] != 3:
        raise ValueError("Image array must have shape (height, width, 3)")

    png_info = PngImagePlugin.PngInfo()

    metadata_dict = _scene_to_metadata_dict(scene)
    for key, value in metadata_dict.items():
        png_info.add_text(key, str(value))

    image_8bit = _dither_to_8bit(image_array)
    image = Image.fromarray(image_8bit, mode="RGB")

    image.save(output_path, "PNG", pnginfo=png_info)


def _dither_to_8bit(image_16bit: np.ndarray) -> np.ndarray:
    """Convert 16-bit RGB to 8-bit with Floyd-Steinberg dithering.

    Args:
        image_16bit: 16-bit RGB array with shape (height, width, 3)

    Returns:
        8-bit RGB array with shape (height, width, 3)
    """
    height, width, channels = image_16bit.shape
    result = image_16bit.copy().astype(np.float32)

    for y in range(height):
        for x in range(width):
            for c in range(channels):
                old_pixel = result[y, x, c]
                new_pixel = np.round(old_pixel / 257.0) * 257.0
                result[y, x, c] = new_pixel

                quant_error = old_pixel - new_pixel

                if x + 1 < width:
                    result[y, x + 1, c] += quant_error * 7.0 / 16.0
                if x - 1 >= 0 and y + 1 < height:
                    result[y + 1, x - 1, c] += quant_error * 3.0 / 16.0
                if y + 1 < height:
                    result[y + 1, x, c] += quant_error * 5.0 / 16.0
                if x + 1 < width and y + 1 < height:
                    result[y + 1, x + 1, c] += quant_error * 1.0 / 16.0

    return (result / 257.0).astype(np.uint8)


def _scene_to_metadata_dict(scene: Scene) -> Dict[str, str]:
    """Convert Scene dataclass to metadata dictionary.

    Args:
        scene: Scene metadata object

    Returns:
        Dictionary with string values for PNG text chunks
    """
    return {
        "body": scene.body.name,
        "utc_time": scene.utc_time,
        "latitude": str(scene.latitude),
        "longitude": str(scene.longitude),
        "altitude_m": str(scene.altitude_m),
        "random_seed": str(scene.random_seed),
        "renderer_id": scene.renderer_id,
    }


def extract_metadata(image_path: str) -> Optional[Dict[str, str]]:
    """Extract metadata from a PNG file.

    Args:
        image_path: Path to PNG file

    Returns:
        Dictionary of metadata key-value pairs, or None if no metadata found
    """
    try:
        image = Image.open(image_path)
    except FileNotFoundError:
        return None

    metadata = {}
    if hasattr(image, "text"):
        for key, value in image.text.items():
            metadata[key] = value

    return metadata if metadata else None
