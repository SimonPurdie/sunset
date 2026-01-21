"""Tests for Software Renderer (SPEC §5.6)."""

import math

import numpy as np
import pytest

from sunset.color.encoding import linear_to_srgb_16bit
from sunset.models.bodies import SOLAR_SYSTEM_BODIES
from sunset.models.observer import Observer
from sunset.renderer import Renderer, render_scene


@pytest.fixture
def earth_observer():
    """Create an Earth observer at sunset."""
    return Observer(
        latitude=40.0,
        longitude=-74.0,
        altitude_m=0.0,
        body=SOLAR_SYSTEM_BODIES["earth"],
        solar_elevation_deg=-0.5,
        sun_direction=(0.0, -0.0087, 0.9999),
        solar_angular_diameter_deg=0.53,
    )


@pytest.fixture
def moon_observer():
    """Create a Moon observer (airless body)."""
    return Observer(
        latitude=0.0,
        longitude=0.0,
        altitude_m=0.0,
        body=SOLAR_SYSTEM_BODIES["moon"],
        solar_elevation_deg=0.0,
        sun_direction=(0.0, 0.0, 1.0),
        solar_angular_diameter_deg=0.52,
    )


class TestRendererBasics:
    """Tests for renderer initialization and basic constraints."""

    def test_minimum_resolution(self):
        """Renderer requires at least 1024 × 512 resolution."""
        with pytest.raises(ValueError, match="Width must be at least 1024"):
            Renderer(width=512, height=256)

        with pytest.raises(ValueError, match="Height must be at least 512"):
            Renderer(width=1024, height=256)

    def test_default_resolution(self):
        """Default resolution is 1024 × 512."""
        renderer = Renderer()
        assert renderer.width == 1024
        assert renderer.height == 512

    def test_custom_resolution(self):
        """Custom resolution must meet minimums."""
        renderer = Renderer(width=2048, height=1024)
        assert renderer.width == 2048
        assert renderer.height == 1024


class TestHorizonGeometry:
    """Tests for horizon presence and geometry (SPEC §5.6 Atomic Tests)."""

    def test_horizon_present(self, earth_observer):
        """Horizon is visible and horizontal."""
        renderer = Renderer()
        image = renderer.render(earth_observer)

        center_col = image.shape[1] // 2
        center_row = image.shape[0] // 2

        above_horizon = image[center_row - 50, center_col]
        below_horizon = image[center_row + 50, center_col]

        below_brightness = np.sum(below_horizon)
        above_brightness = np.sum(above_horizon)

        assert below_brightness < above_brightness * 0.5, (
            "Ground should be darker than sky"
        )

    def test_horizon_horizontal(self, earth_observer):
        """Horizon is horizontal across the image."""
        renderer = Renderer()
        image = renderer.render(earth_observer)

        center_row = image.shape[0] // 2

        left_ground = np.sum(image[center_row + 50, 50])
        right_ground = np.sum(image[center_row + 50, -50])

        left_sky = np.sum(image[center_row - 50, 50])
        right_sky = np.sum(image[center_row - 50, -50])

        ground_similarity = abs(left_ground - right_ground) / (
            left_ground + right_ground
        )
        sky_similarity = abs(left_sky - right_sky) / (left_sky + right_sky)

        assert ground_similarity < 0.5, (
            "Ground brightness should be similar on left and right"
        )
        assert sky_similarity < 0.35, (
            "Sky brightness should be similar on left and right"
        )


class TestSunRendering:
    """Tests for sun disc rendering (SPEC §5.6 Atomic Tests)."""

    def test_sun_intersects_horizon(self, earth_observer):
        """Sun intersects horizon within tolerance."""
        renderer = Renderer()
        image = renderer.render(earth_observer)

        center_col = image.shape[1] // 2
        center_row = image.shape[0] // 2

        sun_pixels = []
        for y in range(center_row - 20, center_row + 20):
            for x in range(center_col - 20, center_col + 20):
                brightness = np.sum(image[y, x])
                if brightness > 50000:
                    sun_pixels.append(y)

        if sun_pixels:
            min_sun_y = min(sun_pixels)
            max_sun_y = max(sun_pixels)
            sun_center_y = (min_sun_y + max_sun_y) / 2

            horizon_row = center_row
            distance_from_horizon = abs(sun_center_y - horizon_row)

            tolerance = image.shape[0] * 0.025

            assert distance_from_horizon < tolerance, (
                f"Sun should be near horizon, but is {distance_from_horizon} pixels away"
            )


class TestNoColorBanding:
    """Tests for color depth and banding (SPEC §5.6 Atomic Tests)."""

    def test_no_color_banding_gradient(self, earth_observer):
        """No visible color banding at 16-bit depth."""
        renderer = Renderer()
        image = renderer.render(earth_observer)

        center_col = image.shape[1] // 2

        vertical_column = image[:, center_col]

        red_values = vertical_column[:, 0]
        green_values = vertical_column[:, 1]
        blue_values = vertical_column[:, 2]

        for channel_values in [red_values, green_values, blue_values]:
            unique_values = np.unique(channel_values)
            ratio = len(unique_values) / len(channel_values)

            assert ratio > 0.7, (
                f"Gradient appears quantized (only {ratio:.2%} unique values)"
            )


class TestAtmosphericBodies:
    """Tests for bodies with atmospheres."""

    def test_earth_sky_colors(self, earth_observer):
        """Earth sunset produces expected warm colors near horizon."""
        renderer = Renderer()
        image = renderer.render(earth_observer)

        center_col = image.shape[1] // 2
        near_horizon = image[image.shape[0] // 2 - 30, center_col]

        near_horizon_sum = np.sum(near_horizon)
        assert near_horizon_sum > 10000, "Near-horizon sky should be bright"

        red_ratio = near_horizon[0] / np.sum(near_horizon)
        green_ratio = near_horizon[1] / np.sum(near_horizon)

        assert red_ratio > green_ratio * 0.8, "Near-horizon should have warm colors"


class TestAirlessBodies:
    """Tests for bodies without atmospheres."""

    def test_moon_black_sky(self, moon_observer):
        """Airless body produces black sky."""
        renderer = Renderer()
        image = renderer.render(moon_observer)

        center_col = image.shape[1] // 2
        high_sky = image[0:50, center_col]

        high_sky_brightness = np.mean(high_sky)
        assert high_sky_brightness < 5000, "High sky should be dark for airless body"


class TestRenderFunction:
    """Tests for convenience render functions."""

    def test_render_scene_function(self, earth_observer):
        """render_scene returns correct shape and type."""
        image = render_scene(earth_observer)

        assert isinstance(image, np.ndarray)
        assert image.shape == (512, 1024, 3)
        assert image.dtype == np.uint16

    def test_render_scene_custom_resolution(self, earth_observer):
        """render_scene respects custom resolution."""
        image = render_scene(earth_observer, width=1280, height=720)

        assert image.shape == (720, 1280, 3)
        assert image.dtype == np.uint16


class TestDeterminism:
    """Tests for deterministic output."""

    def test_deterministic_rendering(self, earth_observer):
        """Same observer produces identical output."""
        renderer = Renderer()

        image1 = renderer.render(earth_observer)
        image2 = renderer.render(earth_observer)

        np.testing.assert_array_equal(image1, image2)
