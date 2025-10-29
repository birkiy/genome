@echo off
REM Minimal Windows build script (may not be used on conda-forge if noarch: python)
python -m pip install . --no-deps -vv
