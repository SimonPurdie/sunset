"""Tests for PNG metadata embedding and extraction."""

import numpy as np
import pytest
from PIL import Image

from sunset.metadata.embedder import embed_metadata, extract_metadata
from sunset.models.bodies import SOLAR_SYSTEM_BODIES
from sunset.models.scene import Scene


@pytest.fixture
def sample_scene():
    """Create a sample scene for testing."""
    return Scene(
        body=SOLAR_SYSTEM_BODIES["earth"],
        utc_time="2025-01-15T18:30:00Z",
        latitude=37.7749,
        longitude=-122.4194,
        altitude_m=0.0,
        random_seed=42,
        renderer_id="test-renderer-1.0.0",
    )


@pytest.fixture
def sample_image_16bit():
    """Create a sample 16-bit RGB image for testing."""
    return np.random.randint(0, 65536, size=(512, 1024, 3), dtype=np.uint16)


@pytest.fixture
def sample_image_8bit():
    """Create a sample 8-bit image for testing (should fail)."""
    return np.random.randint(0, 256, size=(512, 1024, 3), dtype=np.uint8)


def test_metadata_present(tmp_path, sample_scene, sample_image_16bit):
    """Test that all required metadata fields are present in PNG."""
    output_path = tmp_path / "test_output.png"

    embed_metadata(sample_image_16bit, sample_scene, str(output_path))

    metadata = extract_metadata(str(output_path))

    assert metadata is not None, "Metadata should be present"

    required_fields = [
        "body",
        "utc_time",
        "latitude",
        "longitude",
        "altitude_m",
        "random_seed",
        "renderer_id",
    ]

    for field in required_fields:
        assert field in metadata, f"Required field '{field}' missing from metadata"


def test_metadata_parseable(tmp_path, sample_scene, sample_image_16bit):
    """Test that metadata values are parseable to expected types."""
    output_path = tmp_path / "test_output.png"

    embed_metadata(sample_image_16bit, sample_scene, str(output_path))

    metadata = extract_metadata(str(output_path))

    assert metadata is not None

    assert metadata["body"] == "Earth"

    latitude = float(metadata["latitude"])
    assert isinstance(latitude, float)
    assert abs(latitude - 37.7749) < 0.0001

    longitude = float(metadata["longitude"])
    assert isinstance(longitude, float)
    assert abs(longitude - (-122.4194)) < 0.0001

    altitude_m = float(metadata["altitude_m"])
    assert isinstance(altitude_m, float)
    assert abs(altitude_m - 0.0) < 0.0001

    random_seed = int(metadata["random_seed"])
    assert isinstance(random_seed, int)
    assert random_seed == 42

    assert metadata["utc_time"] == "2025-01-15T18:30:00Z"
    assert metadata["renderer_id"] == "test-renderer-1.0.0"


def test_metadata_consistent(tmp_path, sample_scene, sample_image_16bit):
    """Test that metadata matches the rendered scene."""
    output_path = tmp_path / "test_output.png"

    embed_metadata(sample_image_16bit, sample_scene, str(output_path))

    metadata = extract_metadata(str(output_path))

    assert metadata is not None

    assert metadata["body"] == sample_scene.body.name
    assert metadata["utc_time"] == sample_scene.utc_time
    assert float(metadata["latitude"]) == sample_scene.latitude
    assert float(metadata["longitude"]) == sample_scene.longitude
    assert float(metadata["altitude_m"]) == sample_scene.altitude_m
    assert int(metadata["random_seed"]) == sample_scene.random_seed
    assert metadata["renderer_id"] == sample_scene.renderer_id


def test_16bit_png_output(tmp_path, sample_scene, sample_image_16bit):
    """Test that output is 16-bit PNG."""
    output_path = tmp_path / "test_output.png"

    embed_metadata(sample_image_16bit, sample_scene, str(output_path))

    image = Image.open(output_path)
    assert image.mode == "RGB"


def test_invalid_image_type_raises_error(tmp_path, sample_scene, sample_image_8bit):
    """Test that non-uint16 image raises an error."""
    output_path = tmp_path / "test_output.png"

    with pytest.raises(ValueError, match="must be uint16"):
        embed_metadata(sample_image_8bit, sample_scene, str(output_path))


def test_extract_metadata_from_nonexistent_file():
    """Test that extracting metadata from nonexistent file returns None."""
    metadata = extract_metadata("/nonexistent/path/to/file.png")
    assert metadata is None


def test_extract_metadata_from_png_without_metadata(tmp_path, sample_image_16bit):
    """Test extracting metadata from PNG without embedded metadata."""
    output_path = tmp_path / "test_no_metadata.png"

    image_8bit = (sample_image_16bit // 257).astype(np.uint8)
    image = Image.fromarray(image_8bit, mode="RGB")
    image.save(output_path, "PNG")

    metadata = extract_metadata(str(output_path))
    assert metadata is None
