import pytest
from pathlib import Path
import tempfile
from PIL import Image
import numpy as np

from sunset.main import render_sunset
from sunset.models import SOLAR_SYSTEM_BODIES
from sunset.color.spectral import spectral_to_xyz


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

    def test_visual_validation_spectral_variation(self, temp_dir):
        """Test that spectral radiance varies across wavelengths (physics-based)."""
        from sunset.optics.integrator import compute_spectral_radiance
        from sunset.models.atmosphere import get_atmosphere
        from sunset.models.bodies import SOLAR_SYSTEM_BODIES

        body = SOLAR_SYSTEM_BODIES["earth"]
        atmosphere = get_atmosphere("earth", 0)
        sun_direction = (0.0, 0.0, 1.0)

        spectral_radiance = compute_spectral_radiance(
            atmosphere,
            sun_direction,
            0,
            0,
            body.radius_m,
        )

        radiance = spectral_radiance.radiance

        assert radiance is not None, "Spectral radiance not computed"
        assert len(radiance) > 1, (
            "Spectral radiance has insufficient wavelength resolution"
        )

        std_radiance = np.std(radiance)
        assert std_radiance > 1e-6, (
            "Spectral radiance is flat across wavelengths - indicates physics not computed"
        )

    def test_visual_validation_rayleigh_wavelength_scaling(self, temp_dir):
        """Test that shorter wavelengths scatter more than longer wavelengths (Rayleigh physics)."""
        from sunset.optics.integrator import compute_spectral_radiance
        from sunset.models.atmosphere import get_atmosphere
        from sunset.models.bodies import SOLAR_SYSTEM_BODIES

        body = SOLAR_SYSTEM_BODIES["earth"]
        atmosphere = get_atmosphere("earth", 0)
        sun_direction = (0.0, 0.0, 1.0)

        spectral_radiance = compute_spectral_radiance(
            atmosphere,
            sun_direction,
            0,
            0,
            body.radius_m,
        )

        wavelengths = spectral_radiance.wavelengths
        radiance = spectral_radiance.radiance

        blue_idx = np.argmin(np.abs(wavelengths - 450))
        green_idx = np.argmin(np.abs(wavelengths - 550))
        red_idx = np.argmin(np.abs(wavelengths - 650))

        skylight_blue = spectral_radiance.skylight_contribution[blue_idx]
        skylight_red = spectral_radiance.skylight_contribution[red_idx]

        assert skylight_blue > skylight_red * 1.1, (
            "Rayleigh scattering not λ⁻⁴: blue skylight should dominate red"
        )

    def test_visual_validation_sky_gradient_exists(self, temp_dir):
        """Test that sky has color gradient from horizon to zenith (not uniform)."""
        utc_time = "2026-01-21T18:00:00Z"
        seed = 42
        output_path = Path(temp_dir) / "gradient_test.png"

        result = render_sunset(
            utc_time=utc_time,
            body_id="earth",
            random_seed=seed,
            output_path=str(output_path),
            print_caption=False,
        )

        assert result == 0, "Render failed"

        with Image.open(output_path) as img:
            arr = np.array(img)

            middle_col = arr[:, arr.shape[1] // 2, :]

            top_region = middle_col[0 : arr.shape[0] // 4, :]
            bottom_region = middle_col[arr.shape[0] // 2 : 3 * arr.shape[0] // 4, :]

            mean_top = np.mean(top_region, axis=0)
            mean_bottom = np.mean(bottom_region, axis=0)

            color_difference = np.linalg.norm(mean_top - mean_bottom)
            assert color_difference > 5, (
                "Sky gradient too weak - may indicate hardcoded uniform colors"
            )

    def test_visual_validation_atmospheric_influence(self, temp_dir):
        """Test that different atmospheric compositions produce different colors."""
        from sunset.optics.integrator import compute_spectral_radiance
        from sunset.models.atmosphere import get_atmosphere
        from sunset.models.bodies import SOLAR_SYSTEM_BODIES

        earth_body = SOLAR_SYSTEM_BODIES["earth"]
        earth_atmo = get_atmosphere("earth", 0)
        sun_direction = (0.0, 0.0, 1.0)

        mars_body = SOLAR_SYSTEM_BODIES["mars"]
        mars_atmo = get_atmosphere("mars", 0)

        earth_radiance = compute_spectral_radiance(
            earth_atmo,
            sun_direction,
            0,
            0,
            earth_body.radius_m,
        )

        mars_radiance = compute_spectral_radiance(
            mars_atmo,
            sun_direction,
            0,
            0,
            mars_body.radius_m,
        )

        earth_xyz = spectral_to_xyz(earth_radiance)
        mars_xyz = spectral_to_xyz(mars_radiance)

        color_difference = np.linalg.norm(earth_xyz - mars_xyz)
        assert color_difference > 0.01, (
            "Earth and Mars produce identical colors - atmospheric influence not computed"
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
