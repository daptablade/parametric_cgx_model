"""Define imports that are used by one or more modules."""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parents[1]))  # to access the main folder modules
from context import PRECICE_FOLDER
