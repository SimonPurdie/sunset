from dataclasses import dataclass
from typing import Callable, Tuple


@dataclass(frozen=True)
class Body:
    name: str
    radius_m: float
    has_atmosphere: bool
    reference_altitude_m: float = 0.0
    ground_color_base: Tuple[float, float, float] = (0.5, 0.5, 0.5)
    ground_color_variation: float = 0.1


SOLAR_SYSTEM_BODIES = {
    "mercury": Body(
        name="Mercury",
        radius_m=2_439_700,
        has_atmosphere=False,
        ground_color_base=(0.45, 0.43, 0.40),
        ground_color_variation=0.15,
    ),
    "venus": Body(
        name="Venus",
        radius_m=6_051_800,
        has_atmosphere=True,
        ground_color_base=(0.65, 0.55, 0.35),
        ground_color_variation=0.08,
    ),
    "earth": Body(
        name="Earth",
        radius_m=6_371_000,
        has_atmosphere=True,
        ground_color_base=(0.2, 0.3, 0.15),
        ground_color_variation=0.12,
    ),
    "moon": Body(
        name="Moon",
        radius_m=1_737_400,
        has_atmosphere=False,
        ground_color_base=(0.45, 0.45, 0.45),
        ground_color_variation=0.15,
    ),
    "mars": Body(
        name="Mars",
        radius_m=3_389_500,
        has_atmosphere=True,
        ground_color_base=(0.65, 0.35, 0.20),
        ground_color_variation=0.10,
    ),
    "titan": Body(
        name="Titan",
        radius_m=2_574_000,
        has_atmosphere=True,
        ground_color_base=(0.55, 0.50, 0.40),
        ground_color_variation=0.12,
    ),
}
