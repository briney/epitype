"""NanoShaper molecular surface generation wrapper."""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np

from epitype.core.structure import Structure
from epitype.surface.types import SurfacePoint


def generate_surface(
    structure: Structure,
    probe_radius: float = 1.4,
    grid_scale: float = 2.0,
    nanoshaper_path: str | None = None,
) -> list[SurfacePoint]:
    """
    Generate molecular surface using NanoShaper.

    Args:
        structure: Input structure
        probe_radius: Probe radius in Angstroms (default 1.4 for water)
        grid_scale: Grid resolution in grids per Angstrom (default 2.0)
        nanoshaper_path: Path to NanoShaper executable. If None, uses:
            1. EPITYPE_NANOSHAPER_PATH environment variable
            2. Bundled binary (if available for platform)
            3. System PATH lookup

    Returns:
        List of surface points with coordinates and normals
    """
    if nanoshaper_path is None:
        from epitype.bin import get_nanoshaper_path

        nanoshaper_path = get_nanoshaper_path()

    if structure.num_atoms == 0:
        return []

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Write XYZR file (x, y, z, radius for each atom)
        xyzr_path = tmpdir / "structure.xyzr"
        _write_xyzr(structure, xyzr_path)

        # Write NanoShaper parameter file
        prm_path = tmpdir / "config.prm"
        _write_prm_file(prm_path, xyzr_path, probe_radius, grid_scale)

        # Run NanoShaper
        cmd = [nanoshaper_path, str(prm_path)]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=120,
                check=True,
                cwd=str(tmpdir),
            )
        except FileNotFoundError:
            raise RuntimeError(
                f"NanoShaper executable not found at '{nanoshaper_path}'. "
                "Download from: https://gitlab.iit.it/SDecherchi/nanoshaper"
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"NanoShaper failed: {e.stderr}")

        # Parse output files (NanoShaper outputs to working directory)
        vert_path = tmpdir / "triangulatedSurf.vert"
        face_path = tmpdir / "triangulatedSurf.face"

        if not vert_path.exists():
            return []

        surface_points = _parse_nanoshaper_output(vert_path, face_path)

    return surface_points


def _write_xyzr(structure: Structure, path: Path) -> None:
    """Write XYZR file for NanoShaper input."""
    with open(path, "w") as f:
        for atom in structure.atoms:
            f.write(
                f"{atom.coords[0]:.3f} {atom.coords[1]:.3f} "
                f"{atom.coords[2]:.3f} {atom.vdw_radius:.2f}\n"
            )


def _write_prm_file(
    prm_path: Path,
    xyzr_path: Path,
    probe_radius: float,
    grid_scale: float,
) -> None:
    """Write NanoShaper parameter file."""
    content = f"""Grid_scale = {grid_scale}
Grid_perfil = 80.0
Probe_Radius = {probe_radius}
Surface = ses
XYZR_FileName = {xyzr_path}
Build_epsilon_maps = false
Build_status_map = false
Save_Mesh_MSMS_Format = true
Compute_Vertex_Normals = true
"""
    with open(prm_path, "w") as f:
        f.write(content)


def _parse_nanoshaper_output(
    vert_path: Path,
    face_path: Path,
) -> list[SurfacePoint]:
    """Parse NanoShaper vertex and face files (MSMS format)."""
    surface_points: list[SurfacePoint] = []

    # Parse vertex file
    # MSMS format: x y z nx ny nz atom_index ses_area
    with open(vert_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 6:
                continue

            try:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                nx, ny, nz = float(parts[3]), float(parts[4]), float(parts[5])

                # NanoShaper may not include atom_index and area in all modes
                atom_idx = int(parts[6]) if len(parts) > 6 else 0
                area = float(parts[7]) if len(parts) > 7 else 0.0

                point = SurfacePoint(
                    coords=np.array([x, y, z]),
                    normal=np.array([nx, ny, nz]),
                    atom_index=atom_idx,
                    area=area,
                )
                surface_points.append(point)
            except (ValueError, IndexError):
                continue

    return surface_points


def check_nanoshaper_available(nanoshaper_path: str | None = None) -> bool:
    """
    Check if NanoShaper is available.

    Args:
        nanoshaper_path: Path to check. If None, uses default resolution order.
    """
    if nanoshaper_path is None:
        from epitype.bin import get_nanoshaper_path

        nanoshaper_path = get_nanoshaper_path()

    try:
        # NanoShaper prints usage info when run without args and returns non-zero
        subprocess.run(
            [nanoshaper_path],
            capture_output=True,
            timeout=5,
        )
        # NanoShaper returns non-zero when run without args, but that's OK
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False
