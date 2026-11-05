"""Regenerate the exact covariant basis in three steps (requires SymPy)."""
from pathlib import Path
import runpy
ROOT=Path(__file__).resolve().parent
for name in ["build_initial_basis.py","restrict_tangents.py","restrict_family.py"]:
 runpy.run_path(str(ROOT/name),run_name="__main__")
