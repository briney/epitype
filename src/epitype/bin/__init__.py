"""Binary resource management for epitype."""
from __future__ import annotations

import atexit
import os
import platform
import shutil
import stat
import tempfile
from functools import lru_cache
from importlib.resources import as_file, files
from pathlib import Path

# Environment variable for user-specified NanoShaper path
NANOSHAPER_ENV_VAR = "EPITYPE_NANOSHAPER_PATH"

_temp_dir: Path | None = None


def _get_temp_dir() -> Path:
    """Get or create temporary directory for extracted binaries."""
    global _temp_dir
    if _temp_dir is None:
        _temp_dir = Path(tempfile.mkdtemp(prefix="epitype_bin_"))
        atexit.register(lambda: shutil.rmtree(_temp_dir, ignore_errors=True))
    return _temp_dir


def _ensure_executable(path: Path) -> None:
    """Ensure the binary file is executable."""
    current_mode = path.stat().st_mode
    if not (current_mode & stat.S_IXUSR):
        path.chmod(current_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def get_platform_dir() -> str | None:
    """
    Get the platform-specific binary directory name.

    Returns:
        Platform directory name (e.g., 'linux-x86_64') or None if unsupported.
    """
    system = platform.system().lower()
    machine = platform.machine().lower()

    # Map common machine names to standardized names
    machine_map = {
        "x86_64": "x86_64",
        "amd64": "x86_64",
    }
    machine = machine_map.get(machine, machine)

    if system == "linux" and machine == "x86_64":
        return "linux-x86_64"

    return None


@lru_cache(maxsize=1)
def get_bundled_nanoshaper() -> Path | None:
    """
    Get path to bundled NanoShaper binary if available for current platform.

    Extracts the bundled binary to a temporary directory and returns the path.
    Uses caching to ensure single extraction per session.

    Returns:
        Path to bundled binary or None if not available/supported.
    """
    platform_dir = get_platform_dir()
    if platform_dir is None:
        return None

    try:
        bin_package = files("epitype.bin").joinpath(platform_dir)
        with as_file(bin_package.joinpath("NanoShaper")) as src_path:
            if not src_path.exists():
                return None

            # Extract to temp directory for execution
            dest_dir = _get_temp_dir() / platform_dir
            dest_dir.mkdir(parents=True, exist_ok=True)
            dest_path = dest_dir / "NanoShaper"

            if not dest_path.exists():
                shutil.copy2(src_path, dest_path)
                _ensure_executable(dest_path)

            return dest_path
    except (FileNotFoundError, TypeError, OSError):
        return None


def get_nanoshaper_path() -> str:
    """
    Get the NanoShaper path using resolution order.

    Resolution order:
        1. EPITYPE_NANOSHAPER_PATH environment variable
        2. Bundled binary (if available for platform)
        3. System PATH ("NanoShaper")

    Returns:
        Path to NanoShaper executable.
    """
    import warnings

    # 1. Check environment variable
    env_path = os.environ.get(NANOSHAPER_ENV_VAR)
    if env_path:
        path = Path(env_path)
        if path.exists():
            return str(path)
        warnings.warn(
            f"{NANOSHAPER_ENV_VAR} set to '{env_path}' but file not found. "
            "Falling back to bundled or system NanoShaper.",
            stacklevel=2,
        )

    # 2. Check bundled binary
    bundled = get_bundled_nanoshaper()
    if bundled is not None:
        return str(bundled)

    # 3. Fall back to system PATH
    return "NanoShaper"


def get_nanoshaper_info() -> dict[str, str | bool]:
    """
    Get information about NanoShaper configuration.

    Returns:
        Dictionary with keys:
            - path: The resolved path to NanoShaper
            - source: Where the path came from ('environment', 'bundled', 'system')
            - available: Whether NanoShaper is available at the resolved path
            - platform: Current platform string or 'unsupported'
    """
    import subprocess

    env_path = os.environ.get(NANOSHAPER_ENV_VAR)
    bundled = get_bundled_nanoshaper()

    if env_path and Path(env_path).exists():
        path = env_path
        source = "environment"
    elif bundled is not None:
        path = str(bundled)
        source = "bundled"
    else:
        path = "NanoShaper"
        source = "system"

    # Check availability
    try:
        subprocess.run([path], capture_output=True, timeout=5)
        available = True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        available = False

    return {
        "path": path,
        "source": source,
        "available": available,
        "platform": get_platform_dir() or "unsupported",
    }
