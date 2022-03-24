"""Aeroelastic analysis module using preCICE library."""

import os
from pathlib import Path
from multiprocessing import Process
import subprocess

PRECICE_EXECUTES = {
    "aero": "wsl -d Ubuntu-20.04 -e ./wsl_script_aero.bash",
    "structure": "wsl -d Ubuntu-20.04 -e ./wsl_script_structures.bash",
}


def solver_exec(solver, folder):
    """Execute precice call from windows."""

    result = subprocess.run(
        PRECICE_EXECUTES[solver] + "> output_" + solver + ".txt &",
        cwd=folder,
        shell=True,
        check=True,
        capture_output=True,
    )
    print(result)
    print(f"{solver} completed.")


def main():
    """Script that launches the two solvers and records outputs."""

    path_to_precice_modules = Path(os.getcwd(), "precice_modules")
    p1 = Process(target=solver_exec, args=("aero", path_to_precice_modules))
    p1.start()
    p2 = Process(target=solver_exec, args=("structure", path_to_precice_modules))
    p2.start()
    p1.join()
    p2.join()
    print("preCICE run completed")


if __name__ == "__main__":
    main()
