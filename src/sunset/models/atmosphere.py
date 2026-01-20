from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional


@dataclass
class AtmosphericProfile:
    body_name: str
    gas_composition: Dict[str, float]
    rayleigh_coefficient: Callable[[float], float]
    mie_parameters: Optional[Dict[str, float]] = None
    density_profile: Optional[Callable[[float], float]] = None
    pressure_profile: Optional[Callable[[float], float]] = None
    refractive_index: float = 1.0
    absorption_bands: Dict[str, List[tuple[float, float]]] = field(default_factory=dict)


@dataclass
class NullAtmosphere:
    body_name: str


def get_atmosphere(
    body_name: str, altitude_m: float
) -> AtmosphericProfile | NullAtmosphere:
    from ..atmosphere.provider import get_body_profile

    return get_body_profile(body_name, altitude_m)
