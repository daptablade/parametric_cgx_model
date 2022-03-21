import os
from pathlib import Path

# this should match the PRECICE_FOLDER variable defined in precice_modules/aeroelastic_two_way.py
PRECICE_FOLDER = Path(__file__).parent / "precice_modules" / "outputs"
