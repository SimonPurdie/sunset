"""CPU-only software renderer for sunset scenes."""

import math
from typing import Tuple

import numpy as np
from PIL import Image

from ..color.encoding import linear_to_srgb_16bit, xyz_to_srgb
from ..color.spectral import spectral_to_xyz
from ..color.tonemapping import apply_exposure_and_tonemap
from ..models.atmosphere import get_atmosphere
from ..models.observer import Observer
from ..optics.integrator import compute_spectral_radiance


class Renderer:
    """CPU-only sunset renderer."""

    def __init__(
        self, width: int = 1024, height: int = 512, camera_pitch_deg: float = 15.75
    ):
        """Initialize renderer.

        Args:
            width: Image width in pixels (minimum 1024)
            height: Image height in pixels (minimum 512)
            camera_pitch_deg: Camera pitch angle in degrees (positive = looking up, default 15.75° for horizon at 15% from bottom)
        """
        if width < 1024:
            raise ValueError("Width must be at least 1024 pixels")
        if height < 512:
            raise ValueError("Height must be at least 512 pixels")

        self.width = width
        self.height = height
        self.aspect_ratio = width / height
        self.camera_pitch_deg = camera_pitch_deg

        self._precompute_direction_vectors()

    def _precompute_direction_vectors(self):
        """Precompute viewing directions for all pixels (vectorized)."""
        center_x = self.width // 2
        center_y = self.height // 2

        fov_vertical = 90.0
        fov_vertical_rad = math.radians(fov_vertical)
        tan_half_fov = math.tan(fov_vertical_rad / 2.0)

        pitch_rad = math.radians(self.camera_pitch_deg)
        cos_pitch = math.cos(pitch_rad)
        sin_pitch = math.sin(pitch_rad)

        x_coords = np.arange(self.width, dtype=np.float32)
        y_coords = np.arange(self.height, dtype=np.float32)

        X, Y = np.meshgrid(x_coords, y_coords)

        u = (X - center_x) / center_x
        v = (center_y - Y) / center_y

        forward = np.array([0.0, 0.0, 1.0], dtype=np.float32)
        right = np.array([1.0, 0.0, 0.0], dtype=np.float32)
        up = np.array([0.0, 1.0, 0.0], dtype=np.float32)

        rotated_forward = np.array(
            [
                forward[0],
                forward[1] * cos_pitch - forward[2] * sin_pitch,
                forward[1] * sin_pitch + forward[2] * cos_pitch,
            ],
            dtype=np.float32,
        )

        rotated_right = np.array(
            [
                right[0],
                right[1] * cos_pitch - right[2] * sin_pitch,
                right[1] * sin_pitch + right[2] * cos_pitch,
            ],
            dtype=np.float32,
        )

        rotated_up = np.array(
            [
                up[0],
                up[1] * cos_pitch - up[2] * sin_pitch,
                up[1] * sin_pitch + up[2] * cos_pitch,
            ],
            dtype=np.float32,
        )

        dir_x = (
            rotated_forward[0] + u * tan_half_fov * self.aspect_ratio * rotated_right[0]
        )
        dir_y = (
            rotated_forward[1]
            + u * tan_half_fov * self.aspect_ratio * rotated_right[1]
            + v * tan_half_fov * rotated_up[1]
        )
        dir_z = (
            rotated_forward[2]
            + u * tan_half_fov * self.aspect_ratio * rotated_right[2]
            + v * tan_half_fov * rotated_up[2]
        )

        self.direction_map = np.stack([dir_x, dir_y, dir_z], axis=2)

        norms = np.linalg.norm(self.direction_map, axis=2, keepdims=True)
        self.direction_map = self.direction_map / norms

    def render(self, observer: Observer) -> np.ndarray:
        """Render a sunset scene.

        Args:
            observer: Observer location and viewing parameters

        Returns:
            16-bit RGB image array with shape (height, width, 3)
        """
        atmospheric_profile = get_atmosphere(observer.body.name, observer.altitude_m)

        sun_elevation_rad = math.radians(observer.solar_elevation_deg)
        sun_angular_radius = math.radians(observer.solar_angular_diameter_deg / 2.0)

        body_radius_m = observer.body.radius_m
        observer_radius = body_radius_m + observer.altitude_m
        horizon_angle = math.acos(body_radius_m / observer_radius)
        horizon_angle_deg = math.degrees(horizon_angle)

        rgb_image = np.zeros((self.height, self.width, 3), dtype=np.float32)

        v_angles = (
            (self.height // 2 - np.arange(self.height)) / (self.height // 2) * 45.0
        )
        above_horizon_mask = v_angles > -horizon_angle_deg

        sun_direction_enu = np.array(observer.sun_direction, dtype=np.float32)
        sun_direction = np.array(
            [sun_direction_enu[0], sun_direction_enu[1], sun_direction_enu[2]],
            dtype=np.float32,
        )

        dot_products = np.sum(self.direction_map * sun_direction, axis=2)
        dot_products = np.clip(dot_products, -1.0, 1.0)
        angles_to_sun = np.arccos(dot_products)
        sun_mask = angles_to_sun <= sun_angular_radius

        sky_base_srgb = self._compute_sky_base_color(atmospheric_profile, observer)

        sky_brightness_factors = self._compute_sky_brightness_factors(
            self.direction_map, sun_direction
        )

        rgb_image = np.zeros((self.height, self.width, 3), dtype=np.float32)

        for c in range(3):
            sky_channel = (
                sky_base_srgb[c]
                * sky_brightness_factors
                * (0.1 + 0.9 * above_horizon_mask[:, np.newaxis].astype(float))
            )
            rgb_image[:, :, c] = sky_channel

        ground_srgb = self._compute_ground_color(atmospheric_profile, observer)
        for c in range(3):
            rgb_image[~above_horizon_mask, c] = ground_srgb[c]

        sun_color_with_limb_darkening, sun_glow_mask = (
            self._compute_sun_with_limb_darkening_and_glow(
                atmospheric_profile, observer, angles_to_sun, sun_angular_radius
            )
        )
        sun_glow = self._compute_sun_glow(atmospheric_profile, observer, sun_glow_mask)
        rgb_image = np.where(
            sun_mask[:, :, np.newaxis], sun_color_with_limb_darkening, rgb_image
        )
        rgb_image = rgb_image + sun_glow

        rgb_16bit = linear_to_srgb_16bit(rgb_image)

        return rgb_16bit

    def _compute_sky_base_color(
        self, atmospheric_profile, observer: Observer
    ) -> np.ndarray:
        """Compute base sky color (once per scene)."""
        spectral_radiance = compute_spectral_radiance(
            atmospheric_profile,
            observer.sun_direction,
            observer.altitude_m,
            observer.solar_elevation_deg,
            observer.body.radius_m,
        )

        xyz = spectral_to_xyz(spectral_radiance)
        xyz_tonemapped = apply_exposure_and_tonemap(xyz)
        srgb = xyz_to_srgb(xyz_tonemapped)

        return np.clip(srgb, 0.0, 1.0)

    def _compute_sun_color(self, atmospheric_profile, observer: Observer) -> np.ndarray:
        """Compute sun disc color (once per scene)."""
        spectral_radiance = compute_spectral_radiance(
            atmospheric_profile,
            observer.sun_direction,
            observer.altitude_m,
            observer.solar_elevation_deg,
            observer.body.radius_m,
        )

        xyz = spectral_to_xyz(spectral_radiance)
        xyz_tonemapped = apply_exposure_and_tonemap(xyz)
        srgb = xyz_to_srgb(xyz_tonemapped)

        return np.clip(srgb * 1.5, 0.0, 1.0)

    def _compute_ground_color(
        self, atmospheric_profile, observer: Observer
    ) -> np.ndarray:
        """Compute ground color (once per scene)."""
        spectral_radiance = compute_spectral_radiance(
            atmospheric_profile,
            observer.sun_direction,
            observer.altitude_m,
            observer.solar_elevation_deg,
            observer.body.radius_m,
        )

        xyz = spectral_to_xyz(spectral_radiance)
        xyz_tonemapped = apply_exposure_and_tonemap(xyz)
        srgb = xyz_to_srgb(xyz_tonemapped)

        ground_color = srgb * 0.1
        return np.clip(ground_color, 0.0, 1.0)

    def _compute_sun_with_limb_darkening_and_glow(
        self,
        atmospheric_profile,
        observer: Observer,
        angles_to_sun: np.ndarray,
        sun_angular_radius: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute sun disc color with limb darkening effect.

        Limb darkening causes the sun to appear darker at the edges because we're looking
        through cooler outer layers at the limb versus deeper, hotter layers at center.

        Uses simple linear limb darkening law: I(μ) = I₀ * (1 - u * (1 - μ))
        where μ = cos(θ) with θ being angle between line of sight and surface normal.

        Args:
            atmospheric_profile: Atmospheric profile for the body
            observer: Observer parameters
            angles_to_sun: Array of angles from each pixel to sun center
            sun_angular_radius: Angular radius of the sun disc in radians

        Returns:
            Tuple of (sun_color_array with limb darkening, glow_mask for atmospheric glow)
        """
        sun_srgb = self._compute_sun_color(atmospheric_profile, observer)

        mu = np.cos(angles_to_sun)
        limb_darkening_coeff = 0.6

        limb_darkening_factor = np.clip(
            1.0 - limb_darkening_coeff * (1.0 - mu), 0.0, 1.0
        )

        sun_color_with_ld = np.zeros((self.height, self.width, 3), dtype=np.float32)
        for c in range(3):
            sun_color_with_ld[:, :, c] = sun_srgb[c] * limb_darkening_factor

        glow_extent = sun_angular_radius * 2.5
        sun_glow_mask = (angles_to_sun > sun_angular_radius) & (
            angles_to_sun <= glow_extent
        )
        glow_distance_factor = (angles_to_sun - sun_angular_radius) / (
            glow_extent - sun_angular_radius
        )

        return sun_color_with_ld, sun_glow_mask

    def _compute_sun_glow(
        self, atmospheric_profile, observer: Observer, glow_mask: np.ndarray
    ) -> np.ndarray:
        """Compute atmospheric glow/halo effect around the sun disc.

        Args:
            atmospheric_profile: Atmospheric profile for the body
            observer: Observer parameters
            glow_mask: Boolean mask indicating pixels in glow region

        Returns:
            RGB array with glow effect (zeros outside glow region)
        """
        sun_srgb = self._compute_sun_color(atmospheric_profile, observer)

        glow_color = sun_srgb * 0.3
        glow = np.zeros((self.height, self.width, 3), dtype=np.float32)

        angles_to_sun = np.arccos(
            np.clip(
                np.sum(
                    self.direction_map
                    * np.array(observer.sun_direction, dtype=np.float32),
                    axis=2,
                ),
                -1.0,
                1.0,
            )
        )

        sun_angular_radius = math.radians(observer.solar_angular_diameter_deg / 2.0)
        glow_extent = sun_angular_radius * 2.5

        in_glow_region = (angles_to_sun > sun_angular_radius) & (
            angles_to_sun <= glow_extent
        )
        if np.any(in_glow_region):
            normalized_distance = (
                angles_to_sun[in_glow_region] - sun_angular_radius
            ) / (glow_extent - sun_angular_radius)
            fade_factor = 1.0 - normalized_distance
            fade_factor = fade_factor[:, np.newaxis]

            for c in range(3):
                glow[in_glow_region, c] = glow_color[c] * fade_factor[:, 0]

        return glow

    def _compute_sky_brightness_factors(
        self, direction_map: np.ndarray, sun_direction: np.ndarray
    ) -> np.ndarray:
        """Compute brightness factors for sky based on vertical angle from horizon.

        During sunset, sky is brighter near horizon (closer to sun) and darker at zenith.
        The overall brightness scales with solar elevation - lower sun = dimmer overall.
        """
        v_angles = (
            (self.height // 2 - np.arange(self.height)) / (self.height // 2) * 45.0
        )

        brightness_factor_by_height = np.clip(
            0.3 + 0.7 * np.cos(np.radians(v_angles)), 0.1, 1.0
        )

        brightness_factors = np.tile(
            brightness_factor_by_height[:, np.newaxis], (1, self.width)
        )

        return brightness_factors


def render_scene(
    observer: Observer,
    width: int = 1024,
    height: int = 512,
    camera_pitch_deg: float = 15.75,
) -> np.ndarray:
    """Render a sunset scene for the given observer.

    Args:
        observer: Observer location and viewing parameters
        width: Image width in pixels (minimum 1024)
        height: Image height in pixels (minimum 512)
        camera_pitch_deg: Camera pitch angle in degrees (positive = looking up, default 15.75° for horizon at 15% from bottom)

    Returns:
        16-bit RGB image array with shape (height, width, 3)
    """
    renderer = Renderer(width, height, camera_pitch_deg)
    return renderer.render(observer)


def render_scene_to_pil(
    observer: Observer,
    width: int = 1024,
    height: int = 512,
    camera_pitch_deg: float = 15.75,
) -> Image.Image:
    """Render a sunset scene and return as PIL Image.

    Args:
        observer: Observer location and viewing parameters
        width: Image width in pixels (minimum 1024)
        height: Image height in pixels (minimum 512)
        camera_pitch_deg: Camera pitch angle in degrees (positive = looking up, default 15.75° for horizon at 15% from bottom)

    Returns:
        PIL Image in RGB mode
    """
    rgb_array = render_scene(observer, width, height, camera_pitch_deg)

    rgb_8bit = (rgb_array / 257).astype(np.uint8)
    image = Image.fromarray(rgb_8bit, mode="RGB")

    return image
