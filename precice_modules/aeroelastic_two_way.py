from __future__ import division

# ensure packages are available at linux distribution level (not using virtual environment)
from pathlib import Path
import argparse
import numpy as np
import precice
import subprocess
from local_context import PRECICE_FOLDER
from precice_post import write_solver_output_to_file, post_process_solver_iters

SOLVER_PYTHON_PATH = {
    "SolverOne": "venv/Scripts/python.exe",
    "SolverTwo": "venv/Scripts/python.exe",
}
# check that the run folder exist - else create it
if not PRECICE_FOLDER.is_dir():
    PRECICE_FOLDER.mkdir(parents=True)

parser = argparse.ArgumentParser()
parser.add_argument(
    "configurationFileName", help="Name of the xml config file.", type=str
)
parser.add_argument("participantName", help="Name of the solver.", type=str)
parser.add_argument("meshName", help="Name of the mesh.", type=str)

try:
    args = parser.parse_args()
except SystemExit:
    print("")
    print("Usage: python3 aeroelastic_two_way.py config.xml participant-name mesh-name")
    quit()

configuration_file_name = args.configurationFileName
participant_name = args.participantName
mesh_name = args.meshName

if participant_name == "SolverOne":
    write_data_name = "dataOne"
    read_data_name = "dataTwo"

if participant_name == "SolverTwo":
    read_data_name = "dataOne"
    write_data_name = "dataTwo"

solver_process_index = 0
solver_process_size = 1

interface = precice.Interface(
    participant_name, configuration_file_name, solver_process_index, solver_process_size
)

mesh_id = interface.get_mesh_id(mesh_name)
dimensions = interface.get_dimensions()


def read_from_file(file):
    """Read arrays from text files."""
    try:
        array = np.loadtxt(file)
        return array
    except Exception as error:
        print("Could not read solver input file - " + error)


def write_to_file(file, data):
    """Write array to text file."""
    try:
        assert isinstance(
            data, type(np.array([0.0]))
        ), "data should be of type np.ndarray"
        np.savetxt(file, data, fmt="%s")
    except Exception as error:
        print("Could not write solver output file - " + error)


solver_output = []


def execute_python(solver, module):
    """Execute python from local windows venv folder."""
    output = subprocess.run(
        SOLVER_PYTHON_PATH[solver] + " " + module,
        cwd="../.",
        shell=True,
        check=True,
        capture_output=True,
    )
    solver_output.append(output.stdout.decode("utf-8"))


if participant_name == "SolverOne":
    if not (PRECICE_FOLDER / "solver_1_nodes.txt").is_file():
        execute_python(solver=participant_name, module="strip_theory_aerodynamics.py")
    vertices = read_from_file(file=PRECICE_FOLDER / "solver_1_nodes.txt")  # node x,y,z
    read_data = np.zeros((len(vertices), dimensions))  # displacements
    write_data = np.zeros((len(vertices), dimensions))  # forces

if participant_name == "SolverTwo":
    if not (PRECICE_FOLDER / "solver_2_nodes.txt").is_file():
        execute_python(solver=participant_name, module="parametric_box.py")
    vertices = read_from_file(file=PRECICE_FOLDER / "solver_2_nodes.txt")  # node x,y,z
    read_data = np.zeros((len(vertices), dimensions))  # forces
    write_data = np.zeros((len(vertices), dimensions))  # displacements

vertex_ids = interface.set_mesh_vertices(mesh_id, vertices)

# if participant_name == "SolverTwo":
read_data_id = interface.get_data_id(read_data_name, mesh_id)

# if participant_name == "SolverOne":
write_data_id = interface.get_data_id(write_data_name, mesh_id)

dt = interface.initialize()

iter_counter = 0
while interface.is_coupling_ongoing():
    if interface.is_action_required(precice.action_write_iteration_checkpoint()):
        print(f"{participant_name}: Writing iteration checkpoint")
        write_solver_output_to_file(
            solver=participant_name, output=solver_output, folder=PRECICE_FOLDER
        )
        interface.mark_action_fulfilled(precice.action_write_iteration_checkpoint())

    if interface.is_read_data_available():
        read_data = interface.read_block_vector_data(read_data_id, vertex_ids)

    if participant_name == "SolverOne":
        # write displacement input data to file
        write_to_file(file=PRECICE_FOLDER / "solver_1_displacement.txt", data=read_data)
        # execute the aerodynamic analysis with the updated displacements
        execute_python(solver=participant_name, module="strip_theory_aerodynamics.py")
        # update wrtie_data with the force array
        write_data = read_from_file(file=PRECICE_FOLDER / "solver_1_forces.txt")

    if participant_name == "SolverTwo":
        # write force input data to file
        write_to_file(file=PRECICE_FOLDER / "solver_2_forces.txt", data=read_data)
        # execute the aerodynamic analysis with the updated displacements
        execute_python(solver=participant_name, module="parametric_box.py")
        # update wrtie_data with the discplacement array
        write_data = read_from_file(file=PRECICE_FOLDER / "solver_2_displacements.txt")

    if interface.is_write_data_required(dt):
        # if participant_name == "SolverOne":
        interface.write_block_vector_data(write_data_id, vertex_ids, write_data)

    print(f"{participant_name} advancing in time")
    dt = interface.advance(dt)
    iter_counter += 1

    if interface.is_action_required(precice.action_read_iteration_checkpoint()):
        print(f"{participant_name}: Reading iteration checkpoint")
        interface.mark_action_fulfilled(precice.action_read_iteration_checkpoint())

# save iteration history
write_solver_output_to_file(
    solver=participant_name, output=solver_output, folder=PRECICE_FOLDER
)

interface.finalize()
print(f"{participant_name}: Closing python solver ...")
