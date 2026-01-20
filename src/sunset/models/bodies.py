from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class Body:
    name: str
    radius_m: float
    has_atmosphere: bool
    reference_altitude_m: float = 0.0


SOLAR_SYSTEM_BODIES = {
    "mercury": Body(name="Mercury", radius_m=2_439_700, has_atmosphere=False),
    "venus": Body(name="Venus", radius_m=6_051_800, has_atmosphere=True),
    "earth": Body(name="Earth", radius_m=6_371_000, has_atmosphere=True),
    "moon": Body(name="Moon", radius_m=1_737_400, has_atmosphere=False),
    "mars": Body(name="Mars", radius_m=3_389_500, has_atmosphere=True),
    "titan": Body(name="Titan", radius_m=2_574_000, has_atmosphere=True),
}
