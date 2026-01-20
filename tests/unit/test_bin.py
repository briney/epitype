"""Unit tests for binary resource management."""
import os
import stat
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from epitype.bin import (
    NANOSHAPER_ENV_VAR,
    _ensure_executable,
    _get_temp_dir,
    get_bundled_nanoshaper,
    get_nanoshaper_info,
    get_nanoshaper_path,
    get_platform_dir,
)


class TestGetPlatformDir:
    """Tests for get_platform_dir function."""

    def test_returns_linux_x86_64_on_linux_x86(self):
        """Test returns 'linux-x86_64' on Linux x86_64."""
        with patch("epitype.bin.platform.system", return_value="Linux"):
            with patch("epitype.bin.platform.machine", return_value="x86_64"):
                assert get_platform_dir() == "linux-x86_64"

    def test_returns_linux_x86_64_on_amd64(self):
        """Test returns 'linux-x86_64' when machine reports amd64."""
        with patch("epitype.bin.platform.system", return_value="Linux"):
            with patch("epitype.bin.platform.machine", return_value="amd64"):
                assert get_platform_dir() == "linux-x86_64"

    def test_returns_none_on_macos(self):
        """Test returns None on macOS (not yet supported)."""
        with patch("epitype.bin.platform.system", return_value="Darwin"):
            with patch("epitype.bin.platform.machine", return_value="x86_64"):
                assert get_platform_dir() is None

    def test_returns_none_on_windows(self):
        """Test returns None on Windows."""
        with patch("epitype.bin.platform.system", return_value="Windows"):
            with patch("epitype.bin.platform.machine", return_value="AMD64"):
                assert get_platform_dir() is None

    def test_returns_none_on_arm64(self):
        """Test returns None on ARM64 Linux (not yet supported)."""
        with patch("epitype.bin.platform.system", return_value="Linux"):
            with patch("epitype.bin.platform.machine", return_value="aarch64"):
                assert get_platform_dir() is None

    def test_case_insensitive_system(self):
        """Test that system name matching is case-insensitive."""
        with patch("epitype.bin.platform.system", return_value="LINUX"):
            with patch("epitype.bin.platform.machine", return_value="x86_64"):
                assert get_platform_dir() == "linux-x86_64"


class TestGetTempDir:
    """Tests for _get_temp_dir function."""

    def test_creates_directory(self):
        """Test that temp directory is created."""
        import epitype.bin as bin_module

        original_temp_dir = bin_module._temp_dir
        bin_module._temp_dir = None

        try:
            temp_dir = _get_temp_dir()
            assert temp_dir.exists()
            assert temp_dir.is_dir()
            assert "epitype_bin_" in temp_dir.name
        finally:
            bin_module._temp_dir = original_temp_dir

    def test_returns_same_directory_on_repeated_calls(self):
        """Test that repeated calls return the same directory."""
        dir1 = _get_temp_dir()
        dir2 = _get_temp_dir()
        assert dir1 == dir2


class TestEnsureExecutable:
    """Tests for _ensure_executable function."""

    def test_makes_file_executable(self, tmp_path: Path):
        """Test that file is made executable."""
        test_file = tmp_path / "test_binary"
        test_file.write_text("#!/bin/bash\necho test")

        # Remove execute bits
        test_file.chmod(0o644)

        _ensure_executable(test_file)

        mode = test_file.stat().st_mode
        assert mode & stat.S_IXUSR  # User execute
        assert mode & stat.S_IXGRP  # Group execute
        assert mode & stat.S_IXOTH  # Other execute

    def test_preserves_existing_permissions(self, tmp_path: Path):
        """Test that existing permissions are preserved."""
        test_file = tmp_path / "test_binary"
        test_file.write_text("#!/bin/bash\necho test")
        test_file.chmod(0o644)

        _ensure_executable(test_file)

        mode = test_file.stat().st_mode
        assert mode & stat.S_IRUSR  # Read preserved
        assert mode & stat.S_IWUSR  # Write preserved

    def test_noop_if_already_executable(self, tmp_path: Path):
        """Test no change if file is already executable."""
        test_file = tmp_path / "test_binary"
        test_file.write_text("#!/bin/bash\necho test")
        test_file.chmod(0o755)

        original_mode = test_file.stat().st_mode
        _ensure_executable(test_file)

        assert test_file.stat().st_mode == original_mode


class TestGetBundledNanoshaper:
    """Tests for get_bundled_nanoshaper function."""

    def test_returns_none_on_unsupported_platform(self):
        """Test returns None when platform is unsupported."""
        with patch("epitype.bin.get_platform_dir", return_value=None):
            get_bundled_nanoshaper.cache_clear()
            result = get_bundled_nanoshaper()
            assert result is None

    def test_extracts_binary_on_supported_platform(self):
        """Test that binary is extracted on supported platform."""
        get_bundled_nanoshaper.cache_clear()
        result = get_bundled_nanoshaper()

        # Result depends on whether binary exists in package
        # On supported platform with bundled binary, should return path
        if result is not None:
            assert result.exists()
            assert result.name == "NanoShaper"
            # Check it's executable
            assert result.stat().st_mode & stat.S_IXUSR

    def test_caches_result(self):
        """Test that result is cached (lru_cache)."""
        get_bundled_nanoshaper.cache_clear()

        result1 = get_bundled_nanoshaper()
        result2 = get_bundled_nanoshaper()

        # Should be same object (cached)
        assert result1 is result2

    def test_handles_missing_binary_gracefully(self):
        """Test handles case where binary doesn't exist in package."""
        with patch("epitype.bin.get_platform_dir", return_value="nonexistent-platform"):
            get_bundled_nanoshaper.cache_clear()
            result = get_bundled_nanoshaper()
            assert result is None


class TestGetNanoshaperPath:
    """Tests for get_nanoshaper_path function."""

    def test_env_var_takes_priority(self, tmp_path: Path):
        """Test that EPITYPE_NANOSHAPER_PATH env var takes priority."""
        fake_ns = tmp_path / "my_nanoshaper"
        fake_ns.write_text("#!/bin/bash\necho test")

        with patch.dict(os.environ, {NANOSHAPER_ENV_VAR: str(fake_ns)}):
            result = get_nanoshaper_path()
            assert result == str(fake_ns)

    def test_warns_if_env_path_not_found(self, tmp_path: Path):
        """Test warning if env var path doesn't exist."""
        nonexistent = tmp_path / "nonexistent"

        with patch.dict(os.environ, {NANOSHAPER_ENV_VAR: str(nonexistent)}):
            with pytest.warns(UserWarning, match="not found"):
                get_nanoshaper_path()

    def test_falls_back_to_bundled(self, tmp_path: Path):
        """Test falls back to bundled binary if env var not set."""
        # Ensure env var is not set
        env = os.environ.copy()
        env.pop(NANOSHAPER_ENV_VAR, None)

        mock_bundled = tmp_path / "bundled_ns"
        mock_bundled.write_text("#!/bin/bash")

        with patch.dict(os.environ, env, clear=True):
            with patch("epitype.bin.get_bundled_nanoshaper", return_value=mock_bundled):
                result = get_nanoshaper_path()
                assert result == str(mock_bundled)

    def test_falls_back_to_system_path(self):
        """Test falls back to 'NanoShaper' if bundled not available."""
        env = os.environ.copy()
        env.pop(NANOSHAPER_ENV_VAR, None)

        with patch.dict(os.environ, env, clear=True):
            with patch("epitype.bin.get_bundled_nanoshaper", return_value=None):
                result = get_nanoshaper_path()
                assert result == "NanoShaper"

    def test_resolution_order(self, tmp_path: Path):
        """Test the complete resolution order."""
        # Create test files
        env_ns = tmp_path / "env_nanoshaper"
        env_ns.write_text("#!/bin/bash")

        bundled_ns = tmp_path / "bundled_nanoshaper"
        bundled_ns.write_text("#!/bin/bash")

        # Priority 1: Env var
        with patch.dict(os.environ, {NANOSHAPER_ENV_VAR: str(env_ns)}):
            with patch("epitype.bin.get_bundled_nanoshaper", return_value=bundled_ns):
                assert get_nanoshaper_path() == str(env_ns)

        # Priority 2: Bundled (when env var not set)
        env = os.environ.copy()
        env.pop(NANOSHAPER_ENV_VAR, None)
        with patch.dict(os.environ, env, clear=True):
            with patch("epitype.bin.get_bundled_nanoshaper", return_value=bundled_ns):
                assert get_nanoshaper_path() == str(bundled_ns)

        # Priority 3: System (when nothing else available)
        with patch.dict(os.environ, env, clear=True):
            with patch("epitype.bin.get_bundled_nanoshaper", return_value=None):
                assert get_nanoshaper_path() == "NanoShaper"


class TestGetNanoshaperInfo:
    """Tests for get_nanoshaper_info function."""

    def test_returns_dict_with_required_keys(self):
        """Test that info dict contains all required keys."""
        info = get_nanoshaper_info()

        assert "path" in info
        assert "source" in info
        assert "available" in info
        assert "platform" in info

    def test_source_is_environment_when_env_var_set(self, tmp_path: Path):
        """Test source is 'environment' when env var is used."""
        fake_ns = tmp_path / "nanoshaper"
        fake_ns.write_text("#!/bin/bash")

        with patch.dict(os.environ, {NANOSHAPER_ENV_VAR: str(fake_ns)}):
            with patch("subprocess.run"):
                info = get_nanoshaper_info()
                assert info["source"] == "environment"
                assert info["path"] == str(fake_ns)

    def test_source_is_bundled_when_bundled_used(self, tmp_path: Path):
        """Test source is 'bundled' when bundled binary is used."""
        env = os.environ.copy()
        env.pop(NANOSHAPER_ENV_VAR, None)

        mock_bundled = tmp_path / "bundled_ns"
        mock_bundled.write_text("#!/bin/bash")

        with patch.dict(os.environ, env, clear=True):
            with patch("epitype.bin.get_bundled_nanoshaper", return_value=mock_bundled):
                with patch("subprocess.run"):
                    info = get_nanoshaper_info()
                    assert info["source"] == "bundled"

    def test_source_is_system_when_fallback_used(self):
        """Test source is 'system' when falling back to PATH."""
        env = os.environ.copy()
        env.pop(NANOSHAPER_ENV_VAR, None)

        with patch.dict(os.environ, env, clear=True):
            with patch("epitype.bin.get_bundled_nanoshaper", return_value=None):
                info = get_nanoshaper_info()
                assert info["source"] == "system"
                assert info["path"] == "NanoShaper"

    def test_available_is_false_when_not_found(self):
        """Test available is False when binary not found."""
        env = os.environ.copy()
        env.pop(NANOSHAPER_ENV_VAR, None)

        with patch.dict(os.environ, env, clear=True):
            with patch("epitype.bin.get_bundled_nanoshaper", return_value=None):
                with patch("subprocess.run", side_effect=FileNotFoundError):
                    info = get_nanoshaper_info()
                    assert info["available"] is False

    def test_platform_reports_correctly(self):
        """Test platform is reported correctly."""
        with patch("epitype.bin.get_platform_dir", return_value="linux-x86_64"):
            info = get_nanoshaper_info()
            assert info["platform"] == "linux-x86_64"

    def test_platform_reports_unsupported(self):
        """Test platform reports 'unsupported' when not available."""
        with patch("epitype.bin.get_platform_dir", return_value=None):
            info = get_nanoshaper_info()
            assert info["platform"] == "unsupported"
