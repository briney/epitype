"""Geometry utility functions."""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def normalize(v: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Normalize vector to unit length.

    Args:
        v: Input vector

    Returns:
        Unit vector in same direction, or zero vector if input is near-zero
    """
    norm = np.linalg.norm(v)
    if norm < 1e-10:
        return np.zeros_like(v)
    return v / norm


def angle_between(v1: NDArray[np.float64], v2: NDArray[np.float64]) -> float:
    """
    Calculate angle between two vectors in degrees.

    Args:
        v1: First vector
        v2: Second vector

    Returns:
        Angle in degrees (0-180)
    """
    v1_n = normalize(v1)
    v2_n = normalize(v2)

    dot = np.clip(np.dot(v1_n, v2_n), -1.0, 1.0)
    return float(np.degrees(np.arccos(dot)))


def dihedral_angle(
    p1: NDArray[np.float64],
    p2: NDArray[np.float64],
    p3: NDArray[np.float64],
    p4: NDArray[np.float64],
) -> float:
    """
    Calculate dihedral angle between four points in degrees.

    The dihedral angle is the angle between the planes defined by
    (p1, p2, p3) and (p2, p3, p4).

    Args:
        p1: First point
        p2: Second point
        p3: Third point
        p4: Fourth point

    Returns:
        Dihedral angle in degrees (-180 to 180)
    """
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    n1 = normalize(np.cross(b1, b2))
    n2 = normalize(np.cross(b2, b3))

    m1 = np.cross(n1, normalize(b2))

    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    return float(np.degrees(np.arctan2(y, x)))


def rotation_matrix(
    axis: NDArray[np.float64],
    angle: float,
) -> NDArray[np.float64]:
    """
    Create rotation matrix for rotation around axis by angle.

    Uses Rodrigues' rotation formula.

    Args:
        axis: Rotation axis (will be normalized)
        angle: Rotation angle in radians

    Returns:
        3x3 rotation matrix
    """
    axis = normalize(axis)
    c = np.cos(angle)
    s = np.sin(angle)
    t = 1 - c

    x, y, z = axis

    return np.array(
        [
            [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
            [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
            [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
        ]
    )


def centroid(coords: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Calculate centroid of coordinates.

    Args:
        coords: Array of coordinates, shape (N, 3)

    Returns:
        Centroid coordinates, shape (3,)
    """
    return coords.mean(axis=0)


def rmsd(
    coords1: NDArray[np.float64],
    coords2: NDArray[np.float64],
) -> float:
    """
    Calculate RMSD between two coordinate sets.

    Args:
        coords1: First coordinate set, shape (N, 3)
        coords2: Second coordinate set, shape (N, 3)

    Returns:
        Root mean square deviation
    """
    diff = coords1 - coords2
    return float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))
