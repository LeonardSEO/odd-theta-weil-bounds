"""Compatibility entry point for the public read-only certificate checker."""
from pathlib import Path
import runpy

runpy.run_path(str(Path(__file__).resolve().parents[1] / 'scripts/check_certificates.py'), run_name='__main__')
