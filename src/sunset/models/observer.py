from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .bodies import Body


@dataclass(frozen=True)
class Observer:
    latitude: float
    longitude: float
    altitude_m: float
    body: "Body"
    solar_elevation_deg: float
    sun_direction: tuple[float, float, float]
    solar_angular_diameter_deg: float
