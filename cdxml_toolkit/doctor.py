#!/usr/bin/env python3
"""
cdxml-doctor — Diagnostic and test runner for cdxml-toolkit.

Reports dependency availability, runs the test suite, and prints
ChemScript setup instructions when needed.
"""

import os
import sys


def _check_dependency(name, import_fn):
    """Try importing a dependency, return (ok, version_or_error)."""
    try:
        mod = import_fn()
        ver = getattr(mod, "__version__", getattr(mod, "version", "OK"))
        return True, str(ver)
    except Exception as e:
        return False, str(e).split("\n")[0][:60]


def _check_chemdraw_com():
    """Check if ChemDraw COM is registered (without launching it)."""
    try:
        import winreg
        key = winreg.OpenKey(winreg.HKEY_CLASSES_ROOT,
                             "ChemDraw.Application\\CLSID")
        winreg.CloseKey(key)
        return True, "registered"
    except Exception:
        return False, "not registered"


def _check_chemscript():
    """Check if ChemScript bridge is configured and working."""
    try:
        from cdxml_toolkit.chemdraw.chemscript_bridge import (
            ChemScriptBridge, _load_config,
        )
        cfg = _load_config()
        if not cfg.get("python32"):
            return False, "not configured (no python32 in config)"
        if not cfg.get("dll_dir"):
            return False, "not configured (no dll_dir in config)"
        cs = ChemScriptBridge()
        result = cs._call("ping")
        cs.close()
        if result.get("ok"):
            return True, "OK"
        return False, result.get("error", "ping failed")
    except Exception as e:
        return False, str(e).split("\n")[0][:60]


def _print_diagnostics():
    """Print system and dependency diagnostics."""
    print("=== cdxml-toolkit doctor ===\n")

    # System info
    try:
        from cdxml_toolkit import __version__ as pkg_ver
    except (ImportError, AttributeError):
        pkg_ver = "unknown"
    print(f"  Python:       {sys.version.split()[0]} ({sys.executable})")
    print(f"  Platform:     {sys.platform}")
    print(f"  Package:      {pkg_ver}")
    print()

    # Dependencies
    checks = [
        ("RDKit", lambda: __import__("rdkit.Chem", fromlist=["Chem"])),
        ("PyYAML", lambda: __import__("yaml")),
        ("pywin32", lambda: __import__("win32com.client", fromlist=["client"])),
        ("OPSIN", lambda: __import__("py2opsin")),
        ("DECIMER", lambda: __import__("DECIMER")),
        ("pdfplumber", lambda: __import__("pdfplumber")),
        ("MCP", lambda: __import__("mcp")),
        ("pythonnet", lambda: __import__("pythonnet")),
    ]

    print("  Dependencies:")
    for name, fn in checks:
        ok, info = _check_dependency(name, fn)
        status = f"OK ({info})" if ok else f"MISSING ({info})"
        print(f"    {name:14s} {status}")

    # ChemDraw COM
    ok, info = _check_chemdraw_com()
    print(f"    {'ChemDraw COM':14s} {'OK' if ok else 'MISSING'} ({info})")

    # ChemScript
    cs_ok, cs_info = _check_chemscript()
    print(f"    {'ChemScript':14s} {'OK' if cs_ok else 'NOT CONFIGURED'} ({cs_info})")
    print()

    return cs_ok


def _print_chemscript_setup():
    """Print ChemScript setup instructions based on detected DLLs."""
    from cdxml_toolkit.chemdraw.chemscript_bridge import (
        _find_chemdraw_root, _find_chemscript_dlls,
        _find_python_for_chemscript, CONFIG_PATH,
    )

    print("=== ChemScript setup ===\n")

    # Find DLLs
    root = _find_chemdraw_root()
    dll_info = None
    if root:
        dll_info = _find_chemscript_dlls(root)

    if not dll_info:
        print("  ChemScript DLLs not found.")
        print("  ChemDraw with ChemScript must be installed.")
        print("  If installed, run: cdxml-convert --configure")
        return

    bitness = dll_info.get("bitness")
    print(f"  Found ChemScript DLLs:")
    print(f"    Managed:  {dll_info['managed_path']} ({bitness}-bit)")
    print(f"    Native:   {dll_info['native_path']} ({bitness}-bit)")
    print()

    # Check if suitable Python exists
    py = _find_python_for_chemscript(bitness)
    if py:
        print(f"  Python for ChemScript: {py}")
        print(f"  Run: cdxml-convert --configure")
    elif bitness == 32:
        user = os.environ.get("USERNAME", "YOU")
        print("  To enable ChemScript (32-bit DLLs detected):\n")
        print("    set CONDA_SUBDIR=win-32 && conda create -n chemscript32 python=3.10 pip -y")
        print(f"    C:\\Users\\{user}\\miniconda3\\envs\\chemscript32\\python.exe -m pip install pythonnet")
        print("    cdxml-convert --configure")
    elif bitness == 64:
        print("  To enable ChemScript (64-bit DLLs detected):\n")
        print("    pip install pythonnet")
        print("    cdxml-convert --configure")
        print()
        print("  (64-bit ChemScript can use your current Python interpreter)")
    else:
        print("  Could not determine DLL bitness.")
        print("  Run: cdxml-convert --configure")
    print()


def _run_tests():
    """Run the pytest test suite. Returns exit code."""
    print("=== Running tests ===\n")

    # Find tests directory
    try:
        import cdxml_toolkit
        pkg_dir = os.path.dirname(os.path.abspath(cdxml_toolkit.__file__))
        project_root = os.path.dirname(pkg_dir)
        tests_dir = os.path.join(project_root, "tests")
    except Exception:
        tests_dir = None

    if not tests_dir or not os.path.isdir(tests_dir):
        # Try relative to cwd
        tests_dir = os.path.join(os.getcwd(), "tests")

    if not os.path.isdir(tests_dir):
        print("  Tests directory not found.")
        print("  To run tests, cd into the cdxml-toolkit source directory first.")
        print()
        return 1

    try:
        import pytest
        return pytest.main([tests_dir, "-v", "--tb=short"])
    except ImportError:
        print("  pytest not installed. Run: pip install pytest")
        return 1


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(
        description="Diagnostic and test runner for cdxml-toolkit.",
    )
    parser.add_argument("--no-tests", action="store_true",
                        help="Skip running the test suite")
    args = parser.parse_args(argv)

    cs_ok = _print_diagnostics()

    if not args.no_tests:
        exit_code = _run_tests()
        print()
    else:
        exit_code = 0

    if not cs_ok:
        _print_chemscript_setup()

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
