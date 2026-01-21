from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..models.atmosphere import AtmosphericProfile, NullAtmosphere
    from ..models.scene import Scene


def generate_caption(
    scene: "Scene",
    atmosphere: "AtmosphericProfile | NullAtmosphere | None" = None,
) -> str:
    """Generate a caption explaining the sunset scene.

    Args:
        scene: Scene metadata containing body, location, time, etc.
        atmosphere: Optional atmospheric profile for atmospheric explanation.

    Returns:
        Human-readable caption string.
    """
    parts = [f"Sunset on {scene.body.name}"]

    location_str = f"at {abs(scene.latitude):.4f}°{'N' if scene.latitude >= 0 else 'S'}, {abs(scene.longitude):.4f}°{'E' if scene.longitude >= 0 else 'W'}"

    if scene.altitude_m > 0:
        location_str += f", {scene.altitude_m:.0f} meters above surface"

    parts.append(location_str)
    parts.append(f"at {scene.utc_time}")

    if atmosphere is None:
        if not scene.body.has_atmosphere:
            parts.append("No atmosphere")
    else:
        if isinstance(atmosphere, type(None)):
            if not scene.body.has_atmosphere:
                parts.append("No atmosphere")
        else:
            from ..models.atmosphere import NullAtmosphere

            if isinstance(atmosphere, NullAtmosphere):
                parts.append("No atmosphere")
            else:
                primary_gas = max(
                    atmosphere.gas_composition.items(), key=lambda x: x[1]
                )[0]
                parts.append(f"Atmosphere primarily composed of {primary_gas}")

    return ". ".join(parts) + "."
