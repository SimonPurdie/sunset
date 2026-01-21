import math
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..models.atmosphere import AtmosphericProfile, NullAtmosphere
    from ..models.observer import Observer


VISIBILITY_THRESHOLD = 1500.0
MAX_ALTITUDE_M = 200_000.0
ALTITUDE_STEP_M = 1000.0

BODY_SURFACE_PRESSURES = {
    "Venus": 92.0,
    "Earth": 1.0,
    "Mars": 0.006,
    "Titan": 1.5,
}


def calculate_optical_depth(
    profile: "AtmosphericProfile",
    observer_altitude_m: float,
    sun_elevation_deg: float,
) -> float:
    """Calculate optical depth from observer to sun through atmosphere.

    Args:
        profile: Atmospheric profile
        observer_altitude_m: Observer altitude in meters
        sun_elevation_deg: Solar elevation angle in degrees

    Returns:
        Optical depth (dimensionless)
    """
    if profile.density_profile is None:
        return 0.0

    max_sin_angle = 50000.0 / 6_371_000.0
    path_length_factor = 1.0 / math.sin(
        max(math.radians(sun_elevation_deg), max_sin_angle)
    )

    density = profile.density_profile(observer_altitude_m)

    surface_pressure = BODY_SURFACE_PRESSURES.get(profile.body_name, 1.0)

    base_optical_depth = surface_pressure * density * path_length_factor

    if profile.mie_parameters is not None:
        cloud_factor = profile.mie_parameters.get("cloud_concentration", 1.0)
        aerosol_factor = profile.mie_parameters.get("aerosol_concentration", 1.0)
        dust_factor = profile.mie_parameters.get("dust_concentration", 1.0)
        haze_factor = profile.mie_parameters.get("haze_concentration", 1.0)
        mie_factor = cloud_factor + aerosol_factor + dust_factor + haze_factor
        base_optical_depth *= (1.0 + mie_factor) / 2.0

    return base_optical_depth


def resolve_visibility_elevation(
    profile: "AtmosphericProfile | NullAtmosphere",
    observer: "Observer",
) -> tuple[float, str, bool]:
    """Determine if sun is visible and adjust altitude if needed.

    Args:
        profile: Atmospheric profile (could be NullAtmosphere)
        observer: Initial observer location with sun direction

    Returns:
        Tuple of (final_altitude_m, justification, is_visible)
        where is_visible indicates if the body is suitable for rendering
    """
    from ..models.atmosphere import NullAtmosphere

    if isinstance(profile, NullAtmosphere):
        return (
            observer.altitude_m,
            f"{profile.body_name} has no atmosphere; sun always visible",
            True,
        )

    current_altitude = observer.altitude_m
    justification = ""
    is_visible = True

    optical_depth = calculate_optical_depth(
        profile, current_altitude, observer.solar_elevation_deg
    )

    if optical_depth <= VISIBILITY_THRESHOLD:
        justification = (
            f"Sun visible at surface; optical depth {optical_depth:.3f} "
            f"below threshold {VISIBILITY_THRESHOLD}"
        )
        return (current_altitude, justification, is_visible)

    for _ in range(int(MAX_ALTITUDE_M / ALTITUDE_STEP_M)):
        current_altitude += ALTITUDE_STEP_M
        optical_depth = calculate_optical_depth(
            profile, current_altitude, observer.solar_elevation_deg
        )

        if optical_depth <= VISIBILITY_THRESHOLD:
            justification = (
                f"Raised altitude to {current_altitude / 1000:.1f} km; "
                f"optical depth {optical_depth:.3f} below threshold {VISIBILITY_THRESHOLD}"
            )
            return (current_altitude, justification, is_visible)

    is_visible = False
    justification = (
        f"Cannot find visible altitude below {MAX_ALTITUDE_M / 1000:.1f} km; "
        f"optical depth at max altitude {optical_depth:.3f} exceeds threshold {VISIBILITY_THRESHOLD}"
    )

    return (observer.altitude_m, justification, is_visible)
