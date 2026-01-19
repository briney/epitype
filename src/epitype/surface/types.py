"""Surface data types."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass
class SurfacePoint:
    """Single point on molecular surface."""

    coords: NDArray[np.float64]  # (3,) position
    normal: NDArray[np.float64]  # (3,) outward normal
    atom_index: int  # Closest atom
    area: float  # Associated surface area
