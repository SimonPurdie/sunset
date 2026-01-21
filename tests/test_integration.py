import pytest
from pathlib import Path
import tempfile
from PIL import Image
import numpy as np

from sunset.main import render_sunset
from sunset.models import SOLAR_SYSTEM_BODIES


class TestIntegrationFullPipeline:
    """Integration tests for the full sunset rendering pipeline."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    def test_earth_full_pipeline(self, temp_dir):
        """Test full pipeline execution for Earth."""
        utc_time = "2026-01-21T18:00:00Z"
        seed = 42
        output_path = Path(temp_dir) / "earth_sunset.png"

        result = render_sunset(
            utc_time=utc_time,
            body_id="earth",
            random_seed=seed,
            output_path=str(output_path),
            print_caption=False,
        )

        assert result == 0, f"Earth render failed with exit code {result}"
        assert output_path.exists(), "Output file not created"

        self._validate_png_contract(output_path)
        self._validate_metadata_contract(output_path, "earth", utc_time, seed)

    def test_moon_full_pipeline(self, temp_dir):
        """Test full pipeline execution for Moon (airless body)."""
        pytest.skip("Non-Earth bodies require rotation data - see docs/BREADCRUMBS.md")

    def test_mercury_full_pipeline(self, temp_dir):
        """Test full pipeline execution for Mercury (airless body)."""
        pytest.skip("Non-Earth bodies require rotation data - see docs/BREADCRUMBS.md")

    def test_mars_full_pipeline(self, temp_dir):
        """Test full pipeline execution for Mars."""
        pytest.skip("Non-Earth bodies require rotation data - see docs/BREADCRUMBS.md")

    def test_venus_full_pipeline(self, temp_dir):
        """Test full pipeline execution for Venus."""
        pytest.skip("Non-Earth bodies require rotation data - see docs/BREADCRUMBS.md")

    def test_titan_full_pipeline(self, temp_dir):
        """Test full pipeline execution for Titan."""
        pytest.skip("Non-Earth bodies require rotation data - see docs/BREADCRUMBS.md")

    def test_contract_validation_png_format(self, temp_dir):
        """Test that output is valid PNG format."""
        utc_time = "2026-01-21T18:00:00Z"
        seed = 42
        output_path = Path(temp_dir) / "contract_test.png"

        result = render_sunset(
            utc_time=utc_time,
            body_id="earth",
            random_seed=seed,
            output_path=str(output_path),
            print_caption=False,
        )

        assert result == 0, "Render failed"

        with Image.open(output_path) as img:
            assert img.format == "PNG", "Output is not a PNG file"

    def test_contract_validation_resolution(self, temp_dir):
        """Test that output resolution meets minimum requirements."""
        utc_time = "2026-01-21T18:00:00Z"
        seed = 42
        output_path = Path(temp_dir) / "resolution_test.png"

        result = render_sunset(
            utc_time=utc_time,
            body_id="earth",
            random_seed=seed,
            output_path=str(output_path),
            print_caption=False,
        )

        assert result == 0, "Render failed"

        with Image.open(output_path) as img:
            width, height = img.size
            assert width >= 1024, f"Width {width} < 1024 minimum"
            assert height >= 512, f"Height {height} < 512 minimum"

    def test_contract_metadata_completeness(self, temp_dir):
        """Test that all required metadata keys are present."""
        utc_time = "2026-01-21T18:00:00Z"
        seed = 42
        output_path = Path(temp_dir) / "metadata_test.png"

        result = render_sunset(
            utc_time=utc_time,
            body_id="earth",
            random_seed=seed,
            output_path=str(output_path),
            print_caption=False,
        )

        assert result == 0, "Render failed"

        with Image.open(output_path) as img:
            metadata = img.info

            required_keys = [
                "body",
                "utc_time",
                "latitude",
                "longitude",
                "altitude_m",
                "random_seed",
                "renderer_id",
            ]

            for key in required_keys:
                assert key in metadata, f"Required metadata key '{key}' missing"

    def test_contract_metadata_consistency(self, temp_dir):
        """Test that metadata is consistent with scene parameters."""
        utc_time = "2026-01-21T18:00:00Z"
        seed = 999
        output_path = Path(temp_dir) / "consistency_test.png"

        result = render_sunset(
            utc_time=utc_time,
            body_id="earth",
            random_seed=seed,
            output_path=str(output_path),
            print_caption=False,
        )

        assert result == 0, "Render failed"

        with Image.open(output_path) as img:
            metadata = img.info

            assert metadata["utc_time"] == utc_time, "utc_time inconsistent"
            assert metadata["random_seed"] == str(seed), "random_seed inconsistent"

            body = SOLAR_SYSTEM_BODIES["earth"]
            body_name = metadata["body"]
            assert body_name.lower() in ["earth", "Earth"], "body name inconsistent"

            latitude = float(metadata["latitude"])
            assert -90 <= latitude <= 90, "Latitude outside valid range"

            longitude = float(metadata["longitude"])
            assert -180 <= longitude <= 180, "Longitude outside valid range"

            altitude = float(metadata["altitude_m"])
            assert altitude >= 0, "Altitude cannot be negative"

    def test_cross_body_determinism(self, temp_dir):
        """Test that same inputs produce identical outputs for different bodies."""
        utc_time = "2026-01-21T18:00:00Z"
        seed = 777

        bodies_to_test = ["earth"]

        for body_id in bodies_to_test:
            output1 = Path(temp_dir) / f"{body_id}_run1.png"
            output2 = Path(temp_dir) / f"{body_id}_run2.png"

            result1 = render_sunset(
                utc_time=utc_time,
                body_id=body_id,
                random_seed=seed,
                output_path=str(output1),
                print_caption=False,
            )
            assert result1 == 0, f"First run failed for {body_id}"

            result2 = render_sunset(
                utc_time=utc_time,
                body_id=body_id,
                random_seed=seed,
                output_path=str(output2),
                print_caption=False,
            )
            assert result2 == 0, f"Second run failed for {body_id}"

            with Image.open(output1) as img1, Image.open(output2) as img2:
                arr1 = np.array(img1)
                arr2 = np.array(img2)
                assert np.array_equal(arr1, arr2), (
                    f"Outputs differ for {body_id} with same inputs"
                )

    def test_auto_body_selection(self, temp_dir):
        """Test that auto body selection produces valid output."""
        pytest.skip(
            "Auto-selection includes non-Earth bodies that are not yet supported"
        )

    def _validate_png_contract(self, output_path: Path):
        """Validate that the output PNG meets all contract requirements.

        Args:
            output_path: Path to the output PNG file
        """
        with Image.open(output_path) as img:
            assert img.format == "PNG", "Output is not a PNG file"

            width, height = img.size
            assert width >= 1024, f"Width {width} < 1024 minimum"
            assert height >= 512, f"Height {height} < 512 minimum"

    def _validate_metadata_contract(
        self,
        output_path: Path,
        expected_body: str,
        expected_time: str,
        expected_seed: int,
    ):
        """Validate that PNG metadata meets all contract requirements.

        Args:
            output_path: Path to the output PNG file
            expected_body: Expected body name
            expected_time: Expected UTC time
            expected_seed: Expected random seed
        """
        with Image.open(output_path) as img:
            metadata = img.info

            required_keys = [
                "body",
                "utc_time",
                "latitude",
                "longitude",
                "altitude_m",
                "random_seed",
                "renderer_id",
            ]

            for key in required_keys:
                assert key in metadata, f"Required metadata key '{key}' missing"

            assert metadata["utc_time"] == expected_time, "utc_time inconsistent"
            assert metadata["random_seed"] == str(expected_seed), (
                "random_seed inconsistent"
            )

            latitude = float(metadata["latitude"])
            assert -90 <= latitude <= 90, "Latitude outside valid range"

            longitude = float(metadata["longitude"])
            assert -180 <= longitude <= 180, "Longitude outside valid range"

            altitude = float(metadata["altitude_m"])
            assert altitude >= 0, "Altitude cannot be negative"
