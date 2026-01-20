from dataclasses import dataclass, field
from typing import Callable, Dict, Optional


@dataclass
class AtmosphericProfile:
    body_name: str
    gas_composition: Dict[str, float]
    rayleigh_coefficient: Callable[[float], float]
    mie_parameters: Optional[Dict[str, float]] = None
    density_profile: Optional[Callable[[float], float]] = None
    pressure_profile: Optional[Callable[[float], float]] = None
    refractive_index: float = 1.0
    absorption_bands: Dict[str, tuple[float, float]] = field(default_factory=dict)


@dataclass
class NullAtmosphere:
    body_name: str


def get_atmosphere(
    body_name: str, altitude_m: float
) -> AtmosphericProfile | NullAtmosphere:
    from .bodies import SOLAR_SYSTEM_BODIES

    body = SOLAR_SYSTEM_BODIES.get(body_name.lower())
    if body is None:
        raise ValueError(f"Unknown body: {body_name}")

    if not body.has_atmosphere:
        return NullAtmosphere(body_name=body.name)

    raise NotImplementedError(
        f"Atmosphere profiles not yet implemented for {body.name}"
    )
