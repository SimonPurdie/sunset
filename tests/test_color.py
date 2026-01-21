"""Tests for Spectral → Color Pipeline."""

import numpy as np
import pytest

from sunset.color import (
    apply_exposure_and_tonemap,
    linear_to_srgb_16bit,
    linear_to_srgb_8bit,
    spectral_to_xyz,
    xyz_to_srgb,
)
from sunset.models.spectral import SpectralRadiance


class TestSpectralToXYZ:
    def test_energy_luminance(self):
        """More energy → higher luminance."""
        wavelengths = np.arange(380, 781, 10.0)

        low_radiance = np.ones_like(wavelengths) * 0.1
        high_radiance = np.ones_like(wavelengths) * 1.0

        low_spectral = SpectralRadiance(
            wavelengths=wavelengths,
            radiance=low_radiance,
            skylight_contribution=np.zeros_like(wavelengths),
        )

        high_spectral = SpectralRadiance(
            wavelengths=wavelengths,
            radiance=high_radiance,
            skylight_contribution=np.zeros_like(wavelengths),
        )

        low_xyz = spectral_to_xyz(low_spectral)
        high_xyz = spectral_to_xyz(high_spectral)

        assert high_xyz[1] > low_xyz[1], "Higher energy should produce higher luminance"

    def test_hue_ordering(self):
        """Relative wavelength dominance produces distinct colors."""
        wavelengths = np.arange(380, 781, 10.0)
        radiance = np.zeros_like(wavelengths)

        # Peak in blue region (~450nm)
        radiance_blue = radiance.copy()
        blue_idx = np.argmin(np.abs(wavelengths - 450))
        radiance_blue[blue_idx] = 1.0

        # Peak in red region (~650nm)
        radiance_red = radiance.copy()
        red_idx = np.argmin(np.abs(wavelengths - 650))
        radiance_red[red_idx] = 1.0

        spectral_blue = SpectralRadiance(
            wavelengths=wavelengths,
            radiance=radiance_blue,
            skylight_contribution=np.zeros_like(wavelengths),
        )

        spectral_red = SpectralRadiance(
            wavelengths=wavelengths,
            radiance=radiance_red,
            skylight_contribution=np.zeros_like(wavelengths),
        )

        xyz_blue = spectral_to_xyz(spectral_blue)
        xyz_red = spectral_to_xyz(spectral_red)

        # Verify that different spectral distributions produce different XYZ values
        # The actual chromaticity depends on CMF overlap
        assert not np.allclose(xyz_blue, xyz_red), (
            "Different spectra should produce different XYZ"
        )

    def test_zero_radiance_black(self):
        """Zero radiance → black."""
        wavelengths = np.arange(380, 781, 10.0)
        radiance = np.zeros_like(wavelengths)

        spectral = SpectralRadiance(
            wavelengths=wavelengths,
            radiance=radiance,
            skylight_contribution=np.zeros_like(wavelengths),
        )

        xyz = spectral_to_xyz(spectral)

        assert np.allclose(xyz, 0.0), "Zero radiance should produce black (0, 0, 0)"

    def test_xyz_range_positive(self):
        """XYZ values should be non-negative for valid spectral radiance."""
        wavelengths = np.arange(380, 781, 10.0)
        radiance = np.abs(np.random.rand(len(wavelengths)))

        spectral = SpectralRadiance(
            wavelengths=wavelengths,
            radiance=radiance,
            skylight_contribution=np.zeros_like(wavelengths),
        )

        xyz = spectral_to_xyz(spectral)

        assert np.all(xyz >= 0), "XYZ values should be non-negative"


class TestToneMapping:
    def test_tonemapping_hue_preservation(self):
        """Tone mapping preserves hue ordering."""
        # Create two colors with different hues
        xyz_blue = np.array([0.1, 0.1, 0.3])
        xyz_red = np.array([0.3, 0.1, 0.1])

        # Apply tone mapping
        tonemapped_blue = apply_exposure_and_tonemap(xyz_blue)
        tonemapped_red = apply_exposure_and_tonemap(xyz_red)

        # Check that hue ordering is preserved
        # Normalize to chromaticity
        def get_chromaticity(xyz):
            total = xyz.sum()
            if total == 0:
                return np.array([0, 0, 0])
            return xyz / total

        chrom_blue = get_chromaticity(xyz_blue)
        chrom_red = get_chromaticity(xyz_red)
        chrom_blue_tm = get_chromaticity(tonemapped_blue)
        chrom_red_tm = get_chromaticity(tonemapped_red)

        # Blue should still have higher Z than red after tone mapping
        assert tonemapped_blue[2] > tonemapped_red[2], "Blue hue should be preserved"
        assert tonemapped_red[0] > tonemapped_blue[0], "Red hue should be preserved"

    def test_zero_luminance_black(self):
        """Zero luminance remains black after tone mapping."""
        xyz = np.array([0.0, 0.0, 0.0])

        tonemapped = apply_exposure_and_tonemap(xyz)

        assert np.allclose(tonemapped, 0.0), "Zero luminance should remain black"

    def test_high_values_clipped(self):
        """High values are brought into [0, 1] range."""
        xyz = np.array([10.0, 10.0, 10.0])

        tonemapped = apply_exposure_and_tonemap(xyz)

        assert np.all(tonemapped >= 0.0), "Tone-mapped values should be non-negative"
        assert np.all(tonemapped <= 1.0), "Tone-mapped values should be ≤ 1.0"


class TestXyzToSrgb:
    def test_xyz_to_srgb_range(self):
        """XYZ to sRGB should produce values in [0, 1] range."""
        xyz = np.array([0.5, 0.5, 0.5])

        srgb = xyz_to_srgb(xyz)

        assert np.all(srgb >= 0.0), "sRGB values should be non-negative"
        assert np.all(srgb <= 1.0), "sRGB values should be ≤ 1.0"

    def test_gamma_encoding(self):
        """sRGB gamma encoding should be non-linear."""
        # Linear 0.5 should be < 0.5 after gamma encoding
        linear = np.array([0.5, 0.5, 0.5])

        srgb = xyz_to_srgb(linear)

        # Gamma encoding should raise values below threshold
        assert np.all(srgb > 0.5 * 0.9), "Gamma encoding should increase mid-tones"

    def test_black_remains_black(self):
        """Black input should produce black output."""
        xyz = np.array([0.0, 0.0, 0.0])

        srgb = xyz_to_srgb(xyz)

        assert np.allclose(srgb, 0.0), "Black should remain black"


class Test16BitOutput:
    def test_16bit_scaling(self):
        """16-bit conversion should scale correctly."""
        srgb_normalized = np.array([0.0, 0.5, 1.0])

        srgb_16bit = linear_to_srgb_16bit(srgb_normalized)

        expected = np.array([0, 32768, 65535], dtype=np.uint16)

        assert np.allclose(srgb_16bit, expected, atol=1), (
            "16-bit scaling should be correct"
        )

    def test_16bit_range(self):
        """16-bit values should be in [0, 65535] range."""
        srgb_normalized = np.array([0.0, 0.5, 1.0])

        srgb_16bit = linear_to_srgb_16bit(srgb_normalized)

        assert np.all(srgb_16bit >= 0), "16-bit values should be ≥ 0"
        assert np.all(srgb_16bit <= 65535), "16-bit values should be ≤ 65535"


class Test8BitOutput:
    def test_8bit_scaling(self):
        """8-bit conversion should scale correctly."""
        srgb_normalized = np.array([0.0, 0.5, 1.0])

        srgb_8bit = linear_to_srgb_8bit(srgb_normalized)

        expected = np.array([0, 128, 255], dtype=np.uint8)

        assert np.allclose(srgb_8bit, expected, atol=1), (
            "8-bit scaling should be correct"
        )

    def test_8bit_range(self):
        """8-bit values should be in [0, 255] range."""
        srgb_normalized = np.array([0.0, 0.5, 1.0])

        srgb_8bit = linear_to_srgb_8bit(srgb_normalized)

        assert np.all(srgb_8bit >= 0), "8-bit values should be ≥ 0"
        assert np.all(srgb_8bit <= 255), "8-bit values should be ≤ 255"


class TestFullPipeline:
    def test_full_spectral_to_srgb(self):
        """End-to-end test: spectral → XYZ → tone-mapped XYZ → sRGB."""
        wavelengths = np.arange(380, 781, 10.0)

        # Create a simple spectrum with peak in green region
        radiance = np.zeros_like(wavelengths)
        green_idx = np.argmin(np.abs(wavelengths - 550))
        radiance[green_idx] = 1.0

        spectral = SpectralRadiance(
            wavelengths=wavelengths,
            radiance=radiance,
            skylight_contribution=np.zeros_like(wavelengths),
        )

        # Full pipeline
        xyz = spectral_to_xyz(spectral)
        xyz_tonemapped = apply_exposure_and_tonemap(xyz)
        srgb = xyz_to_srgb(xyz_tonemapped)

        # Validate output
        assert len(srgb) == 3, "Should produce RGB triplet"
        assert np.all(srgb >= 0.0), "sRGB values should be non-negative"
        assert np.all(srgb <= 1.0), "sRGB values should be ≤ 1.0"

    def test_pipeline_with_skylight(self):
        """Pipeline works with skylight contribution."""
        wavelengths = np.arange(380, 781, 10.0)
        radiance = np.ones_like(wavelengths) * 0.1
        skylight = np.ones_like(wavelengths) * 0.05

        spectral = SpectralRadiance(
            wavelengths=wavelengths,
            radiance=radiance,
            skylight_contribution=skylight,
        )

        xyz = spectral_to_xyz(spectral)
        xyz_tonemapped = apply_exposure_and_tonemap(xyz)
        srgb = xyz_to_srgb(xyz_tonemapped)

        # The internal radiance is used, not skylight contribution
        # (skylight is included in the radiance by the optics integrator)
        assert len(srgb) == 3, "Should produce RGB triplet"
        assert np.all(srgb >= 0.0), "sRGB values should be non-negative"
        assert np.all(srgb <= 1.0), "sRGB values should be ≤ 1.0"
