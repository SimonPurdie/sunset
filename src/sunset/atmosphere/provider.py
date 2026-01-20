from typing import Literal

from ..models.atmosphere import AtmosphericProfile, NullAtmosphere, get_atmosphere
from ..models.bodies import SOLAR_SYSTEM_BODIES
from .profiles import earth_profile, mars_profile, venus_profile, titan_profile


def get_body_profile(
    body_name: str, altitude_m: float = 0.0
) -> AtmosphericProfile | NullAtmosphere:
    """Get atmospheric profile for a Solar System body.

    Args:
        body_name: Name of the body (case-insensitive)
        altitude_m: Altitude in meters above reference surface

    Returns:
        AtmosphericProfile if body has atmosphere, NullAtmosphere otherwise

    Raises:
        ValueError: If body_name is not recognized
    """
    body_key = body_name.lower()
    body = SOLAR_SYSTEM_BODIES.get(body_key)

    if body is None:
        raise ValueError(f"Unknown body: {body_name}")

    if not body.has_atmosphere:
        return NullAtmosphere(body_name=body.name)

    body_profiles: dict[str, AtmosphericProfile] = {
        "earth": earth_profile,
        "mars": mars_profile,
        "venus": venus_profile,
        "titan": titan_profile,
    }

    profile = body_profiles.get(body_key)

    if profile is None:
        raise ValueError(f"Atmosphere profile not yet implemented for {body.name}")

    return profile
