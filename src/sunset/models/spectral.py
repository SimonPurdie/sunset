from dataclasses import dataclass
import numpy as np


@dataclass
class SpectralRadiance:
    wavelengths: np.ndarray
    radiance: np.ndarray
    skylight_contribution: np.ndarray

    def __post_init__(self):
        if len(self.wavelengths) != len(self.radiance):
            raise ValueError("wavelengths and radiance must have same length")
        if len(self.wavelengths) != len(self.skylight_contribution):
            raise ValueError(
                "wavelengths and skylight_contribution must have same length"
            )
