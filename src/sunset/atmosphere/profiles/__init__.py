from dataclasses import dataclass

from .earth import earth_profile
from .mars import mars_profile
from .venus import venus_profile
from .titan import titan_profile

__all__ = [
    "earth_profile",
    "mars_profile",
    "venus_profile",
    "titan_profile",
]
