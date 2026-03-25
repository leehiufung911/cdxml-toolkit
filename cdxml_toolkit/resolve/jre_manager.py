"""Manage a bundled JRE for OPSIN.

A JRE zip is shipped inside the package (``cdxml_toolkit/_jre/``).
On first use it is extracted to ``~/.cdxml-toolkit/jre/`` so that
py2opsin can run without requiring the user to install Java.

Discovery order (used by :func:`get_java`):
  1. System Java on PATH or JAVA_HOME
  2. Already-extracted JRE at ``~/.cdxml-toolkit/jre/``
  3. Extract from bundled zip (one-time, ~45 MB)
  4. Download from Adoptium API (network fallback)
"""

from __future__ import annotations

import io
import os
import shutil
import sys
import zipfile
from pathlib import Path
from typing import Optional

# Where the extracted JRE lives
_JRE_BASE = Path.home() / ".cdxml-toolkit" / "jre"

# Bundled JRE zip inside the package
_BUNDLED_ZIP = Path(__file__).resolve().parent.parent / "_jre" / "temurin-21-jre-win-x64.zip"

# Network fallback URL (Adoptium API)
_ADOPTIUM_URL = (
    "https://api.adoptium.net/v3/binary/latest/21/ga/windows/x64/jre/"
    "hotspot/normal/eclipse?project=jdk"
)

# Cached result
_java_exe: Optional[str] = None


def _find_system_java() -> Optional[str]:
    """Check PATH and JAVA_HOME for an existing java executable."""
    java = shutil.which("java")
    if java:
        return java

    java_home = os.environ.get("JAVA_HOME")
    if java_home:
        for name in ("java.exe", "java"):
            candidate = os.path.join(java_home, "bin", name)
            if os.path.isfile(candidate):
                return candidate
    return None


def _find_extracted_java() -> Optional[str]:
    """Check ~/.cdxml-toolkit/jre/ for an already-extracted JRE."""
    if not _JRE_BASE.is_dir():
        return None
    # The zip extracts to a subdirectory like jdk-21.0.10+7-jre/
    for entry in _JRE_BASE.iterdir():
        if entry.is_dir():
            for name in ("java.exe", "java"):
                candidate = entry / "bin" / name
                if candidate.is_file():
                    return str(candidate)
    return None


def _extract_bundled_jre() -> Optional[str]:
    """Extract the JRE zip shipped inside the package.

    Returns the path to java.exe, or None if the bundled zip is missing.
    """
    if not _BUNDLED_ZIP.is_file():
        return None

    _JRE_BASE.mkdir(parents=True, exist_ok=True)

    print("  [cdxml-toolkit] Extracting bundled JRE (one-time)...",
          file=sys.stderr)
    try:
        with zipfile.ZipFile(_BUNDLED_ZIP) as zf:
            zf.extractall(_JRE_BASE)
    except Exception as e:
        print(f"  [cdxml-toolkit] JRE extraction failed: {e}", file=sys.stderr)
        return None

    java = _find_extracted_java()
    if java:
        print(f"  [cdxml-toolkit] JRE ready: {java}", file=sys.stderr)
    return java


def _download_jre() -> Optional[str]:
    """Download Eclipse Temurin JRE 21 from Adoptium (network fallback).

    Only used if the bundled zip is missing (e.g. minimal source install).
    Returns the path to java.exe, or None on failure.
    """
    try:
        import urllib.request
    except ImportError:
        return None

    _JRE_BASE.mkdir(parents=True, exist_ok=True)

    print("  [cdxml-toolkit] Downloading JRE for OPSIN (~45 MB)...",
          file=sys.stderr)
    try:
        req = urllib.request.Request(
            _ADOPTIUM_URL,
            headers={"User-Agent": "cdxml-toolkit/0.5"},
        )
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = resp.read()
    except Exception as e:
        print(f"  [cdxml-toolkit] JRE download failed: {e}", file=sys.stderr)
        return None

    print("  [cdxml-toolkit] Extracting JRE...", file=sys.stderr)
    try:
        with zipfile.ZipFile(io.BytesIO(data)) as zf:
            zf.extractall(_JRE_BASE)
    except Exception as e:
        print(f"  [cdxml-toolkit] JRE extraction failed: {e}", file=sys.stderr)
        return None

    java = _find_extracted_java()
    if java:
        print(f"  [cdxml-toolkit] JRE ready: {java}", file=sys.stderr)
    return java


def get_java(download: bool = True) -> Optional[str]:
    """Return the path to a ``java`` executable.

    Discovery order:
      1. System Java (PATH / JAVA_HOME)
      2. Already-extracted JRE at ``~/.cdxml-toolkit/jre/``
      3. Extract from bundled zip (ships with the package)
      4. Download from Adoptium API (network fallback)

    Args:
        download: If True (default), allow network download as last
                  resort when the bundled zip is also missing.

    Returns:
        Absolute path to ``java`` or ``java.exe``, or None.
    """
    global _java_exe
    if _java_exe is not None:
        return _java_exe

    # 1. System Java
    _java_exe = _find_system_java()
    if _java_exe:
        return _java_exe

    # 2. Already-extracted bundled JRE
    _java_exe = _find_extracted_java()
    if _java_exe:
        return _java_exe

    # 3. Extract from bundled zip
    _java_exe = _extract_bundled_jre()
    if _java_exe:
        return _java_exe

    # 4. Network fallback
    if download:
        _java_exe = _download_jre()
        return _java_exe

    return None


def ensure_java_on_path(download: bool = True) -> bool:
    """Make sure ``java`` is discoverable by subprocess calls.

    Finds (or extracts/downloads) a JRE, then adds its ``bin/``
    directory to ``PATH`` and sets ``JAVA_HOME`` so that py2opsin's
    ``subprocess.run(["java", ...])`` works.

    Returns True if Java is available, False otherwise.
    """
    java = get_java(download=download)
    if not java:
        return False

    java_bin_dir = os.path.dirname(java)
    path = os.environ.get("PATH", "")
    if java_bin_dir not in path:
        os.environ["PATH"] = java_bin_dir + os.pathsep + path
    os.environ["JAVA_HOME"] = os.path.dirname(java_bin_dir)
    return True
