from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .bodies import Body


@dataclass(frozen=True)
class Scene:
    body: "Body"
    utc_time: str
    latitude: float
    longitude: float
    altitude_m: float
    random_seed: int
    renderer_id: str
