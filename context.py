"""Define imports that are used by one or more modules."""

from pathlib import Path

# folder accessed by precice and all coupled solvers
PRECICE_FOLDER = Path(__file__).parent / "precice_modules" / "outputs"
