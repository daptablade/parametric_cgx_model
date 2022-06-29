"""
Parametric shell half-wing model geometry creation and meshing in CGX,
followed by FEA analysis (using NASTRAN or CALCULIX).

The script assumes that CGX, CCX and nastran are installed and working.
Adjust the executable calls as needed in LOCAL_EXECUTES.
Nastran or CalculiX CrunchiX are only required if you want to run FEM analyses
(see References in README file).

Execute the python script 'parametric_box.py' and inspect outputs:
choose between
>> main(INPUT[0]) for a metallic Nastran model
>> main(INPUT[1]) for a composite Calculix model
>> main(INPUT[2]) for a multisection composite Calculix model with core
"""

# import external libraries
import os
import shutil
from pathlib import Path
from math import ceil
import warnings
import csv
import subprocess
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation

from context import PRECICE_FOLDER

# set these execution paths to local binaries as available
LOCAL_EXECUTES = {
    "CGX": "wsl cgx_2.19",
    "CALCULIX": "wsl ccx_2.19",
    "NASTRAN": "nastran",  # set to None if not installed
}

PLOT_FLAG = False

# parameters from up-stream processes and analyses
INPUTS = [
    {
        "span": 2.0,
        "chord": 0.2,
        "filled_sections_flags": False,
        "inputs_folder": "inputs",
        "output_folder": Path(os.getcwd(), "outputs"),
        "airfoil_csv_file": "naca0012.csv",
        "analysis_file": "normal_modes.bdf",
        "nele_foil": [10, 10],
        "nele_span": 40,
        "node_merge_tol": 0.002,
        "cgx_ele_type": 9,  # 9: S4, 10: S8 (linear or quadratic elements)
        "cgx_solver": "nas",  # or "nas"
        "fea_solver": "NASTRAN",  # or "NASTRAN"
        "mesh_file": "all.bdf",
        "boundary_conditions": {"fix_lines": [0, 1], "loaded_lines": None},
        "process_flags": {"run_post": False, "delete_run_folder": False},
    },
    {
        "span": 2.0,
        "chord": 0.2,
        "filled_sections_flags": False,
        "inputs_folder": "inputs",
        "output_folder": Path(os.getcwd(), "outputs"),
        "airfoil_csv_file": "naca0012.csv",
        "analysis_file": "ccx_static_tip_shear",  # specify without file extension for CALCULIX
        "nele_foil": [10, 10],
        "nele_span": 40,
        "node_merge_tol": 0.002,
        "cgx_ele_type": 10,  # 9: S4, 10: S8 (linear or quadratic elements)
        "cgx_solver": "abq",  # or "nas"
        "fea_solver": "CALCULIX",  # or "NASTRAN"
        "composite_plies": [
            {
                "id": "p_0",
                "thickness": 0.0002,
                "material": "EL",
                "orientation": "ORI_0",
            },
            {
                "id": "p_90",
                "thickness": 0.0002,
                "material": "EL",
                "orientation": "ORI_90",
            },
        ],
        "orientations": [
            {"id": "ORI_0", "1": [0.0, 1.0, 0.0], "2": [-1.0, 0.0, 0.0]},
            {"id": "ORI_90", "1": [1.0, 0.0, 0.0], "2": [0.0, 1.0, 0.0]},
        ],
        "composite_layup": {
            "aero": (["p_90"] + ["p_0"] * 3 + ["p_90"] * 2 + ["p_0"] * 3 + ["p_90"]),
        },
        "shell_set_name": {"aero": "Eall"},
        "composite_props_file": "composite_shell.inp",
        "mesh_file": "all.msh",
        "boundary_conditions": {"fix_lines": [0, 1], "loaded_lines": [5, 6]},
        "process_flags": {"run_post": True, "delete_run_folder": False},
    },
    {
        "span": [0.085, 0.507, 0.085, 0.646, 0.085, 0.507, 0.085],
        "chord": [0.38, 0.38, 0.38, 0.38, 0.38, 0.38, 0.38],
        "filled_sections_flags": [True, False, True, False, True, False, True],
        "inputs_folder": "inputs",
        "output_folder": Path(os.getcwd(), "outputs"),
        "airfoil_csv_file": "naca4418.csv",
        "airfoil_cut_chord_percentages": [5, 95],
        "analysis_file": "ccx_normal_modes",  # specify without file extension for CALCULIX
        "nele_foil": [4, 10, 2, 2, 10, 4],
        "nele_span": [4, 10, 4, 10, 4, 10, 4],
        "node_merge_tol": 0.00002,
        "cgx_ele_type": 10,  # 9: S4, 10: S8 (linear or quadratic elements)
        "cgx_solver": "abq",  # or "nas"
        "fea_solver": "CALCULIX",  # or "NASTRAN"
        "composite_plies": [
            {
                "id": "glass_0",
                "thickness": 0.0005,
                "material": "Glass_ply",
                "orientation": "ORI_0",
            },
            {
                "id": "dummy",
                "thickness": 0.00025,
                "material": "Foam",
                "orientation": "ORI_0",
            },
        ],
        "orientations": [{"id": "ORI_0", "1": [0.0, 1.0, 0.0], "2": [-1.0, 0.0, 0.0]}],
        "composite_layup": {"ribs": ["dummy"], "aero": ["glass_0"]},
        "shell_set_name": {"ribs": "ERIBS", "aero": "EAERO"},
        "composite_props_file": "composite_shell.inp",
        "mesh_file": "all.msh",
        "boundary_conditions": {"fix_lines": None, "loaded_lines": None},
        "process_flags": {"run_post": False, "delete_run_folder": False},
    },
    {
        "span": 2.0,
        "chord": 0.2,
        "filled_sections_flags": False,
        "inputs_folder": "inputs",
        "output_folder": Path(os.getcwd(), "outputs"),
        "precice_folder": PRECICE_FOLDER,
        "airfoil_csv_file": "naca0012.csv",
        "analysis_file": "ccx_static_aero_forces",  # specify without file extension for CALCULIX
        "nele_foil": [10, 10],
        "nele_span": 40,
        "node_merge_tol": 0.002,
        "cgx_ele_type": 10,  # 9: S4, 10: S8 (linear or quadratic elements)
        "cgx_solver": "abq",  # or "nas"
        "fea_solver": "CALCULIX",  # or "NASTRAN"
        "composite_plies": [
            {
                "id": "p_0",
                "thickness": 0.0002,
                "material": "EL",
                "orientation": "ORI_0",
            },
            {
                "id": "p_90",
                "thickness": 0.0002,
                "material": "EL",
                "orientation": "ORI_90",
            },
        ],
        "orientations": [
            {"id": "ORI_0", "1": [0.0, 1.0, 0.0], "2": [-1.0, 0.0, 0.0]},
            {"id": "ORI_90", "1": [1.0, 0.0, 0.0], "2": [0.0, 1.0, 0.0]},
        ],
        "composite_layup": {
            "aero": (["p_90"] + ["p_0"] * 3 + ["p_90"] * 2 + ["p_0"] * 3 + ["p_90"]),
        },
        "shell_set_name": {"aero": "Eall"},
        "composite_props_file": "composite_shell.inp",
        "mesh_file": "all.msh",
        "boundary_conditions": {
            "fix_lines": [0, 1],
            "loaded_lines": [5],
            "loaded_surfaces": [0],
        },
        "process_flags": {
            "run_post": True,
            "delete_run_folder": False,
            "precice_inout": True,
        },
    },
]


def main(inputs):
    """Create a box FEM model from parametric inputs."""

    print("The parametric inputs are:")
    print(inputs)

    # create the FEM analysis folder and copy input files into it
    run_folder = _make_analysis_folder(
        inputs=[inputs["analysis_file"]],
        inputs_folder=inputs["inputs_folder"],
        outputs_folder=inputs["output_folder"],
    )

    geometry = get_geometry(inputs, plot_flag=PLOT_FLAG)
    infile = get_cgx_input_file(geometry, inputs, run_folder)
    print(f"Created cgx input file {infile}.")

    execute_cgx(infile)
    print("Created analysis input files with CGX.")

    if "composite_layup" in inputs:
        get_composite_properties_input(inputs, run_folder)

    if (
        "precice_inout" in inputs["process_flags"]
        and inputs["process_flags"]["precice_inout"]
    ):
        # check that the run folder exist - else create it
        if not inputs["precice_folder"].is_dir():
            inputs["precice_folder"].mkdir(parents=True)

        if (
            PLOT_FLAG
            or not Path(inputs["precice_folder"], "solver_2_nodes.txt").is_file()
            or not Path(inputs["precice_folder"], "solver_2_node_ids.txt").is_file()
        ):
            node_ids, nodes_xyz = _write_nodes_file(
                readf={
                    "node_ids": Path(run_folder, "TOP.nam"),
                    "nodes": Path(run_folder, inputs["mesh_file"]),
                },
                writef={
                    "nodes": Path(inputs["precice_folder"], "solver_2_nodes.txt"),
                    "node_ids": Path(inputs["precice_folder"], "solver_2_node_ids.txt"),
                },
            )
            if not Path(inputs["precice_folder"], "solver_2_forces.txt").is_file():
                _write_zero_forces_dummy(
                    len(node_ids),
                    writef=Path(inputs["precice_folder"], "solver_2_forces.txt"),
                )

        forces = _write_force_input(
            readf={
                "forces": Path(inputs["precice_folder"], "solver_2_forces.txt"),
                "node_ids": Path(inputs["precice_folder"], "solver_2_node_ids.txt"),
            },
            writef=Path(run_folder, "AERO_FORCES.inp"),
        )
        if PLOT_FLAG:
            # plot force distribution over the nodes
            _plot_forces_ditribution(nodes_xyz, forces)

    # run the FEM model analysis
    execute_fea(inputs["analysis_file"], inputs["fea_solver"], run_folder)
    print("Executed FEM analysis.")

    if inputs["process_flags"]["run_post"]:
        # recover the analysis results
        outputs = get_fea_outputs(
            file=inputs["analysis_file"],
            solver=inputs["fea_solver"],
            mesh_file=inputs["mesh_file"],
            folder=run_folder,
        )
        print(outputs)

    if (
        "precice_inout" in inputs["process_flags"]
        and inputs["process_flags"]["precice_inout"]
    ):
        _write_displacement_output(
            readf=run_folder / (inputs["analysis_file"] + ".dat"),
            writef=Path(inputs["precice_folder"], "solver_2_displacements.txt"),
        )

    if inputs["process_flags"]["delete_run_folder"]:
        # delete the run folder (recommend setting to True for optimisation)
        shutil.rmtree(run_folder)
        print(f"deleted run folder: {run_folder}")

    print("End main process.\n")

    if "outputs" in locals():
        return outputs


def get_geometry(inputs, plot_flag=False):
    """Translate parameters into geometry description that CGX can understand,
    that's points, lines and surfaces."""

    if "airfoil_cut_chord_percentages" not in inputs:
        inputs["airfoil_cut_chord_percentages"] = None

    aerofoil = _get_aerofoil_from_file(
        Path(inputs["inputs_folder"], inputs["airfoil_csv_file"]),
        plot_flag=plot_flag,
        splitpc=inputs["airfoil_cut_chord_percentages"],
    )
    points, seqa, split_points = _get_cgx_points_3d(
        aerofoil, inputs["chord"], inputs["span"]
    )
    lines, rib_surfaces, aero_surfaces, bodies, aero_surfaces_flip = _get_cgx_lines_3d(
        seqa,
        nele_foil=inputs["nele_foil"],
        nele_span=inputs["nele_span"],
        split_points=split_points,
        filled_sections=inputs["filled_sections_flags"],
    )

    return {
        "aerofoil": aerofoil,
        "points": points,
        "point_seqa": seqa,
        "lines": lines,
        "surfaces": {
            "ribs": rib_surfaces,
            "aero": aero_surfaces,
            "aero_surfaces_flip": aero_surfaces_flip,
        },
        "bodies": bodies,
    }


def get_cgx_input_file(geometry, inputs, folder):
    """Write CGX batch commands to file."""

    fdb_geom_file = folder / "cgx_infile.fdb"

    if "boundary_conditions" in inputs:
        fix_lines = inputs["boundary_conditions"]["fix_lines"]  # [0] is to fix the root
        loaded_lines = inputs["boundary_conditions"]["loaded_lines"]
        if "loaded_surfaces" in inputs["boundary_conditions"]:
            loaded_surfaces = inputs["boundary_conditions"]["loaded_surfaces"]
        else:
            loaded_surfaces = None

    else:
        fix_lines = None
        loaded_lines = None
        loaded_surfaces = None

    # create string of all input commands
    cgx_commands = _get_commands(
        geometry,
        fix_lines,
        loaded_lines,
        loaded_surfaces=loaded_surfaces,
        merge_tol=inputs["node_merge_tol"],
        cgx_ele_type=inputs["cgx_ele_type"],
        solver=inputs["cgx_solver"],
    )

    # write string of commands to file
    with open(fdb_geom_file, "w", encoding="utf-8") as f:
        f.write("".join(cgx_commands))

    return fdb_geom_file


def get_composite_properties_input(inputs, run_folder):
    """write an FEA input file with the composite properties."""

    if inputs["fea_solver"] == "CALCULIX":

        # check and update the element type in the mesh input file
        str_find = "*ELEMENT, TYPE=S8,"
        str_replace = "*ELEMENT, TYPE=S8R,"
        _file_find_replace(
            file=(run_folder / inputs["mesh_file"]),
            find=str_find,
            replace_with=str_replace,
        )

        if "filled_sections_flags" in inputs and not isinstance(
            inputs["filled_sections_flags"], list
        ):
            inputs["filled_sections_flags"] = [inputs["filled_sections_flags"]]

        shell_set_name = inputs["shell_set_name"]
        if "filled_sections_flags" in inputs and any(inputs["filled_sections_flags"]):

            assert (
                isinstance(inputs["airfoil_cut_chord_percentages"], list)
                and len(inputs["airfoil_cut_chord_percentages"]) == 2
            ), (
                "if 'filled_sections_flags' is switched on, 'airfoil_cut_chord_percentages'"
                "should be a list of length 2."
            )

            # create separate element sets for shells and solids
            str_find = "*ELEMENT, TYPE=S8R, ELSET=Eall"
            str_replace = "*ELEMENT, TYPE=S8R, ELSET=SURF"
            _file_find_replace(
                file=(run_folder / inputs["mesh_file"]),
                find=str_find,
                replace_with=str_replace,
            )
            str_find = "*ELEMENT, TYPE=C3D20, ELSET=Eall"
            str_replace = "*ELEMENT, TYPE=C3D20, ELSET=CORE"
            _file_find_replace(
                file=(run_folder / inputs["mesh_file"]),
                find=str_find,
                replace_with=str_replace,
            )

        # get input file cards for this solver
        ccx_commands = _get_ccx_composite_shell_props(
            plies=inputs["composite_plies"],
            orientations=inputs["orientations"],
            layup=inputs["composite_layup"],
            shell_set_name=shell_set_name,
        )
    else:
        warnings.warn(
            f"Composite properties not implemented for solver option {inputs['fea_solver']}."
        )

    # write string of commands to file
    with open(run_folder / inputs["composite_props_file"], "w", encoding="utf-8") as f:
        f.write("".join(ccx_commands))


def execute_cgx(infile):
    """Run CGX with the batch input file to generate the mesh output files."""
    if LOCAL_EXECUTES["CGX"]:
        subprocess.run(
            LOCAL_EXECUTES["CGX"] + " -bg " + infile.parts[-1],
            cwd=infile.parent,
            shell=True,
            check=True,
            capture_output=True,
        )
    else:
        raise ValueError("Need to specify an execution path for CalculiX GraphiX.")


def execute_fea(file, solver, run_folder):
    """Run CGX with the batch input file to generate the mesh output files."""

    if LOCAL_EXECUTES[solver]:
        subprocess.run(
            LOCAL_EXECUTES[solver] + " " + file,
            cwd=run_folder,
            shell=True,
            check=True,
            capture_output=True,
        )
    else:
        warnings.warn(
            f"{solver} execution path has not been set. Model analysis is skipped."
        )


def get_fea_outputs(file, solver, mesh_file, folder):
    """Recover the analysis outputs and process them for plotting."""

    if solver == "CALCULIX":
        # read file and recover node displacements
        output_file = folder / (str(file) + ".dat")
        # all_disp : [nodeid, vx, vy, vz]
        all_disps = _get_from_dat(output_file, data="displacements")

        # get node ids for displacement output > read from LAST.nam file
        node_ids = np.array(
            _get_from_inp(file=Path(folder, "LAST.nam"), keyword="NSET,NSET="),
            dtype=float,
        )

        # filter all_disp by node ids
        disp_filter = np.isin(all_disps[:, 0], node_ids)
        tip_disps = all_disps[disp_filter, :]

        # average the deflections
        v_mean = np.zeros((3))
        for disp in range(3):
            v_mean[disp] = np.average(tip_disps[:, disp + 1])

        # calculate the average rotations
        rotations_mean = _get_average_rotation(
            all_disps=tip_disps, mesh_file=folder / mesh_file
        )

        outputs = {
            "Ux": v_mean[0],
            "Uy": v_mean[1],
            "Uz": v_mean[2],
            "Rx": rotations_mean[0],
            "Ry": rotations_mean[1],
            "Rz": rotations_mean[2],
        }
        return outputs

    else:
        warnings.warn(f"Output processing not implemented for solver option {solver}.")
        return {}


########### Private functions that do not get called directly


def _get_aerofoil_from_file(file, plot_flag=True, splitpc=None, pt_offset=6):
    """
    This function reads an aerofoil geometry from csv and calculates tc_max.

    Args:
        file: file with xy-positions of the airfoil outline in Selig format.
              For example http://airfoiltools.com/airfoil/seligdatfile?airfoil=n0012-il

    Returns:
        airfoil data
    """

    # read aerofoil input file
    airfoil = []
    with open(file, mode="r", encoding="utf-8") as infile:
        reader = csv.reader(infile, skipinitialspace=True)
        for row in reader:
            airfoil.append(row)
    name = airfoil[0]
    coordinates = np.array([string[0].split() for string in airfoil[1:]], dtype=float)

    # replace the last coordinate to close the airfoil at the trailing-edge
    coordinates[-1] = coordinates[0]

    # we assume that there is a [0.0, 0.0] point in the airfoil
    LE_index = np.where(coordinates[:, 0] == 0.0)[0][0]
    leading_edge_pt = LE_index

    splits = []
    if splitpc:
        # check that there are enough points to split the section
        min_points = 100
        if len(coordinates) < min_points:
            raise ValueError(
                "The parameter 'airfoil_cut_chord_percentages' requires "
                f"at least {min_points:d} airfoil spline points in 'airfoil_csv_file'"
            )

        # re-order the pc from TE to LE
        splitpc.sort(reverse=True)

        # trim points that are within min number of points form leading or trailing edge
        trimmed_coords = np.hstack(
            [np.array([np.arange(len(coordinates))]).T, coordinates]
        )
        trimmed_coords = np.vstack(
            [
                trimmed_coords[pt_offset : int(LE_index - ceil(pt_offset / 2)), :],
                trimmed_coords[int(LE_index + ceil(pt_offset / 2)) : -pt_offset, :],
            ]
        )
        # find two points that match the percentage chord closely
        for split_number, split in enumerate(splitpc):
            point_distances_x = np.abs(trimmed_coords[:, 1] - split / 100)
            pt = {"top": 0, "bot": 0}
            dist_top = point_distances_x[0]
            dist_bot = point_distances_x[-1]
            for index, dist in enumerate(point_distances_x):
                if dist < dist_top and trimmed_coords[index, 2] > 0:
                    pt["top"] = int(trimmed_coords[index, 0])
                    dist_top = dist
                if dist < dist_bot and trimmed_coords[index, 2] < 0:
                    pt["bot"] = int(trimmed_coords[index, 0])
                    dist_bot = dist

            if split_number >= 1:
                # check number of points separating splits
                if (
                    np.abs(pt["top"] - splits[-1]["top"]) < pt_offset
                    or np.abs(pt["bot"] - splits[-1]["bot"]) < pt_offset
                ):
                    raise ValueError(
                        f"Values {splitpc[split_number-1]} and {split:f} in "
                        "'airfoil_cut_chord_percentages' are too close together."
                    )
            splits.append(pt)

    if plot_flag:
        plt.plot(coordinates[:, 0], coordinates[:, 1], "-xr")
        if splitpc:
            for split in splits:
                plt.plot(
                    [coordinates[split["top"], 0], coordinates[split["bot"], 0]],
                    [coordinates[split["top"], 1], coordinates[split["bot"], 1]],
                    "-b",
                )
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(name)
        plt.show()

    return dict(
        name=name,
        coordinates=coordinates,
        splits=splits,
        leading_edge_pt=leading_edge_pt,
    )


def _get_cgx_points_3d(aerofoil, chord, span):
    """This function generates the CGX input file points and point sequences."""

    if not isinstance(span, list):
        span = [span]
    if not isinstance(chord, list):
        chord = [chord]

    def wing_with_splits(aerofoil, chord, span):
        seqa = []
        starting_y = 0
        pt_counter = 0
        split_points = np.empty((len(aerofoil["splits"]), 2, 0), dtype=int)
        for section_index, length_y in enumerate(span):

            x = aerofoil["coordinates"][:, 0] * chord[section_index]
            z = aerofoil["coordinates"][:, 1] * chord[section_index]

            if section_index == 0:  # only needed at the root of the wing
                y_root = np.ones(x.size) * starting_y
                points = np.vstack([x, y_root, z]).T

            # section tip
            y_tip = np.ones(x.size) * (starting_y + length_y)
            points = np.append(points, np.vstack([x, y_tip, z]).T, axis=0)

            def airfoil_seqa(pt_counter, seqa, all_split_points, x_size):
                # SEQA for first airfoil top splines
                pt = 0
                split_points = []
                for split in aerofoil["splits"]:
                    indices = np.arange(pt_counter + pt + 1, pt_counter + split["top"])
                    seqa.append(indices)
                    pt = split["top"]
                    split_points.append(pt_counter + split["top"])
                # SEQA for first airfoil LE spline
                seqa.append(
                    np.arange(
                        pt_counter + pt + 1,
                        pt_counter + aerofoil["leading_edge_pt"],
                    )
                )

                if any(aerofoil["splits"]):
                    seqa.append(
                        np.arange(
                            pt_counter + aerofoil["leading_edge_pt"] + 1,
                            pt_counter + aerofoil["splits"][-1]["bot"],
                        )
                    )
                    # SEQA for first airfoil bot spline
                    pt = x_size - 1
                    bot_seqa = []
                    for split in aerofoil["splits"]:
                        indices = np.flipud(
                            np.arange(
                                pt_counter + pt - 1, pt_counter + split["bot"], -1
                            )
                        )
                        bot_seqa.append(indices)
                        pt = split["bot"]
                        split_points.append(pt_counter + split["bot"])
                    bot_seqa.reverse()
                    seqa += bot_seqa
                    # all_split_point is nested list of dim 3: aerofoil -> split -> point
                    all_split_points = np.dstack(
                        [
                            all_split_points,
                            np.reshape(split_points, (len(aerofoil["splits"]), 2)).T,
                        ]
                    )
                else:
                    seqa.append(
                        np.arange(
                            pt_counter + aerofoil["leading_edge_pt"] + 1,
                            pt_counter + x_size - 2,
                        )
                    )
                    all_split_points = []

                return seqa, all_split_points

            seqa, split_points = airfoil_seqa(
                pt_counter=pt_counter,
                seqa=seqa,
                all_split_points=split_points,
                x_size=x.size,
            )
            if section_index == 0:  # only needed at the root of the wing
                seqa, split_points = airfoil_seqa(
                    pt_counter=pt_counter + x.size,
                    seqa=seqa,
                    all_split_points=split_points,
                    x_size=x.size,
                )

            starting_y += length_y
            pt_counter = seqa[-1][-1] + 2

        return points, seqa, split_points

    points, seqa, split_points = wing_with_splits(aerofoil, chord, span)

    return points, seqa, split_points


def _get_cgx_lines_3d(
    seqa,
    nele_foil=20,
    nele_span=40,
    nele_split=4,
    split_points=None,
    filled_sections=None,
):
    """This function creates the aerofoil section splines and the spanwise bounding lines in CGX."""

    nele_multiplier = 2  # =2 to account for quadratic elements
    lines = []
    aero_surfaces = []
    aero_surfaces_flip = []
    rib_surfaces = []

    if not isinstance(nele_span, list):
        nele_span = [nele_span]

    if not isinstance(nele_foil, list):
        nele_foil = [nele_foil]

    if not isinstance(filled_sections, list):
        filled_sections = [filled_sections]

    splits = 0
    seqas_per_aerofoil = 2
    if isinstance(split_points, np.ndarray):
        splits = split_points.shape[0]
        seqas_per_aerofoil = splits * 2 + 2
    aerofoils = int(len(seqa) / seqas_per_aerofoil)

    airfoil_index = 0
    lcounter = 0
    for seqa_id, seq in enumerate(seqa):

        # aerofoil lines
        lines.append(
            [
                seq[0] - 1,
                seq[-1] + 1,
                seqa_id,
                nele_foil[seqa_id % seqas_per_aerofoil] * nele_multiplier,
            ]
        )
        lcounter += 1

        if (seqa_id + 1) % seqas_per_aerofoil == 0:

            if isinstance(split_points, np.ndarray):
                for split_index, split in enumerate(split_points[:, :, airfoil_index]):

                    # aerofoil split lines
                    lines.append([split[0], split[1], nele_split * nele_multiplier])
                    lcounter += 1

                    if split_index > 0:
                        # prepare rib surfaces definition
                        rib_surfaces.append(
                            [
                                lcounter - 1 - seqas_per_aerofoil,
                                lcounter - 1,
                                lcounter - split_index * 2 - 2,
                                -(lcounter - 2),
                            ]
                        )

            # spanwise lines at trailing edge
            if (seqa_id + 1) / seqas_per_aerofoil < aerofoils:

                for te_line_inc in range(seqas_per_aerofoil + 1):
                    if te_line_inc < seqas_per_aerofoil:
                        start_id = seqa_id + 1 - seqas_per_aerofoil + te_line_inc
                        end_id = seqa_id + 1 + te_line_inc
                        side = 0
                        pt_offset = -1
                    else:
                        start_id = seqa_id - seqas_per_aerofoil + te_line_inc
                        end_id = seqa_id + te_line_inc
                        side = -1
                        pt_offset = 1
                    lines.append(
                        [
                            seqa[start_id][side] + pt_offset,
                            seqa[end_id][side] + pt_offset,
                            nele_span[airfoil_index] * nele_multiplier,
                        ]
                    )
                    lcounter += 1

                    if te_line_inc < seqas_per_aerofoil:
                        if te_line_inc < seqas_per_aerofoil / 2:  # top surface
                            aero_surfaces_flip.append(True)
                        else:  # bot surface
                            aero_surfaces_flip.append(False)
                        # prepare aero surfaces definition
                        aero_surfaces.append(
                            [
                                lcounter - 1 - splits - seqas_per_aerofoil,
                                lcounter,
                                -(lcounter + seqas_per_aerofoil),
                                -(lcounter - 1),
                            ]
                        )

            airfoil_index += 1

    # check that ptB_id > ptA_id
    assert all(
        [line[0] < line[1] for line in lines]
    ), "something has gone wrong in the line definition."

    # solid bodies
    bodies = []
    for surf_id, _ in enumerate(rib_surfaces[1:]):
        if filled_sections[surf_id]:
            bodies.append([surf_id, surf_id + 1])

    return lines, rib_surfaces, aero_surfaces, bodies, aero_surfaces_flip


def _get_commands(
    geometry,
    fix_lines,
    loaded_lines,
    loaded_surfaces,
    merge_tol=0.001,
    cgx_ele_type=10,
    solver="abq",
    max_entries_per_line=9,
):
    def divide_chunks(l, n):
        # looping till length l
        for i in range(0, len(l), n):
            yield l[i : i + n]

    commands = []

    # points
    for entity_id, point in enumerate(geometry["points"]):
        commands.append(
            f"PNT P{entity_id:05d} {point[0]:e} {point[1]:e} {point[2]:e}\n"
        )

    commands.append("# =============== \n")
    # point sequences
    for entity_id, points in enumerate(geometry["point_seqa"]):
        commands.append(f"SEQA A{entity_id:05d} pnt ")
        for ii in range(0, len(points), 8):
            line_end = " = \n" if ii + 8 < len(points) else "\n"
            commands.append(
                " ".join([f"P{point:05d}" for point in points[ii : ii + 8]]) + line_end
            )

    commands.append("# =============== \n")
    # lines
    for entity_id, line in enumerate(geometry["lines"]):
        if len(line) == 3:  # straight line
            commands.append(
                f"LINE L{entity_id:05d} P{line[0]:05d} P{line[1]:05d} {line[2]:d} \n"
            )
        elif len(line) == 4:  # spline
            commands.append(
                f"LINE L{entity_id:05d} P{line[0]:05d} P{line[1]:05d} A{line[2]:05d} {line[3]:d} \n"
            )

    commands.append("# =============== \n")
    # surfaces
    rib_ids = []
    for entity_id, surf in enumerate(geometry["surfaces"]["ribs"]):
        commands.append(
            f"GSUR V{entity_id:05d} + BLEND "
            + " ".join(
                [
                    f"+ L{np.abs(line):05d}"
                    if np.sign(line) >= 0
                    else f"- L{np.abs(line):05d}"
                    for line in surf
                ]
            )
            + "\n"
        )
        rib_ids.append(entity_id)

    aero_ids = []
    flip_surfaces = []
    for counter, surf in enumerate(geometry["surfaces"]["aero"]):
        entity_id = counter + (rib_ids[-1] if rib_ids else -1) + 1
        commands.append(
            f"GSUR V{entity_id:05d} + BLEND "
            + " ".join(
                [
                    f"+ L{np.abs(line):05d}"
                    if np.sign(line) >= 0
                    else f"- L{np.abs(line):05d}"
                    for line in surf
                ]
            )
            + "\n"
        )
        if geometry["surfaces"]["aero_surfaces_flip"][counter]:
            flip_surfaces.append(f"FLIP V{entity_id:05d}" + "\n")
        aero_ids.append(entity_id)

    commands.append("# =============== \n")
    # bodies
    for entity_id, body in enumerate(geometry["bodies"]):
        commands.append(f"BODY B{entity_id:05d} V{body[0]:05d} V{body[1]:05d}" + "\n")

    commands.append("# =============== \n")
    # SPC and load sets
    if fix_lines:
        for chunk in divide_chunks(fix_lines, max_entries_per_line):
            commands.append(
                "SETA SPC l " + " ".join([f"L{line:05d}" for line in chunk]) + "\n"
            )
    if loaded_lines:
        for chunk in divide_chunks(loaded_lines, max_entries_per_line):
            commands.append(
                "SETA LAST l " + " ".join([f"L{line:05d}" for line in chunk]) + "\n"
            )
    if loaded_surfaces:
        for chunk in divide_chunks(loaded_surfaces, max_entries_per_line):
            commands.append(
                "SETA TOP s " + " ".join([f"V{id:05d}" for id in chunk]) + "\n"
            )

    commands.append("# =============== \n")
    # surface meshes
    surfaces = geometry["surfaces"]["ribs"] + geometry["surfaces"]["aero"]
    for entity_id, _ in enumerate(surfaces):
        commands.append(f"MSHP V{entity_id:05d} s {cgx_ele_type:d} 0 1.000000e+00\n")

    commands.append("")
    # sets of surfaces
    if rib_ids:
        for chunk in divide_chunks(rib_ids, max_entries_per_line):
            commands.append(
                "SETA RIBS s " + " ".join([f"V{id:05d}" for id in chunk]) + "\n"
            )
    if aero_ids:
        for chunk in divide_chunks(aero_ids, max_entries_per_line):
            commands.append(
                "SETA AERO s " + " ".join([f"V{id:05d}" for id in chunk]) + "\n"
            )

    commands.append("# =============== \n")
    # body meshes
    if geometry["bodies"]:
        for entity_id, _ in enumerate(geometry["bodies"]):
            commands.append(f"MSHP B{entity_id:05d} b 4 0 1.000000e+00\n")

    commands.append("# =============== \n")
    # custom export statement
    commands.append("mesh all\n")
    commands.append(f"merg n all {merge_tol:6f} 'nolock'\n")
    commands.append("comp nodes d\n")
    if flip_surfaces:
        commands += flip_surfaces
    if fix_lines:
        commands.append("comp SPC d\n")
        commands.append(f"send SPC {solver} spc 123456\n")
    if loaded_lines:
        commands.append("comp LAST d\n")
        commands.append(f"send LAST {solver} names\n")
    if loaded_surfaces:
        commands.append("comp TOP d\n")
        commands.append(f"send TOP {solver} names\n")
    if rib_ids:
        commands.append("comp RIBS d\n")
        commands.append(f"send RIBS {solver} names\n")
    if aero_ids:
        commands.append("comp AERO d\n")
        commands.append(f"send AERO {solver} names\n")
    commands.append(f"send all {solver} \n")
    commands.append("quit\n")

    return commands


def _get_ccx_composite_shell_props(
    plies=None, orientations=None, layup=None, shell_set_name=None
):

    commands = []
    if not shell_set_name:
        shell_set_name = {"ribs": "ERIBS", "aero": "EAERO"}

    # orientation cards
    for ori in orientations:
        commands.append(f"*ORIENTATION,NAME={ori['id']}\n")
        commands.append(", ".join(str(x) for x in [*ori["1"], *ori["2"]]) + "\n")

    commands.append("** =============== \n")
    # shell property
    for (key, section_name) in shell_set_name.items():
        commands.append(f"*SHELL SECTION,ELSET={section_name},COMPOSITE\n")
        for ply in layup[key]:
            props = [p for p in plies if p["id"] == ply][0]
            commands.append(
                f"{props['thickness']:6f},,{props['material']},{props['orientation']}\n"
            )

    return commands


def _file_find_replace(file, find: str, replace_with: str):
    with open(file, "r", encoding="utf-8") as f:
        contents = f.readlines()

    for index, line in enumerate(contents):
        if find in line:
            contents[index] = line.replace(find, replace_with)
            print(f"Find & Replace edited file '{file}' at line {index:d}.")
            break

    with open(file, "w", encoding="utf-8") as f:
        f.write("".join(contents))


def _get_from_dat(file, data: str = None):
    with open(file, "r", encoding="utf-8") as f:
        contents = f.readlines()

    # extract data assuming this is the only output
    for index, line in enumerate(contents):
        if data in line:
            data_values = contents[index + 2 :]
            break

    # parse data
    if data == "displacements":
        output = np.array(
            [
                [float(val) for val in line.strip().split()]
                for line in data_values
                if line
            ],
            dtype=float,
        )

    return output


def _get_average_rotation(
    all_disps: np.ndarray = None,
    mesh_file: str = None,
    axes: tuple = (
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ),  # basic x, y, z axes
):
    # for each output node, recover the node undeformed mesh position
    # node_positions : [nodeid, x, y, z]
    node_positions = _get_nodes_from_inp(mesh_file, ref_nodes=all_disps[:, 0])

    origin_index = np.where(node_positions[:, 1] == 0.0)[0][0]  # LE
    # calculate the approximate rotation from each node
    rotations = _get_rot_from_disp(all_disps, node_positions, axes, origin_index)

    # average the rotation values
    rotations_mean = np.zeros((3))
    for rot in rotations:
        rotations_mean[rot] = np.average(rotations[rot])

    return rotations_mean


def _get_nodes_from_inp(file, ref_nodes=None):
    with open(file, "r", encoding="utf-8") as f:
        contents = f.readlines()

    # find node definitions in the file
    start = None
    nodes = []
    for index, line in enumerate(contents):
        if "*NODE" in line:
            start = index + 1
        elif start and "*" in line:  # marks the next keyword
            break
        elif start and index >= start:
            nodes.append(line)

    # parse to array
    nodes_vals = np.array(
        [[float(val) for val in line.strip().split(",")] for line in nodes if line],
        dtype=float,
    )

    # filter reference nodes
    if any(ref_nodes):
        nodes_vals = np.array(
            [line for line in nodes_vals if any(line[0] == ref_nodes)], dtype=float
        )

    return nodes_vals


def _get_rot_from_disp(all_disps, node_positions, axes, origin_index):
    """Return an array of rotations in radians. Using small deflection assumptions."""

    all_rots = {key: [] for key in range(len(axes))}
    origin_displacements = all_disps[origin_index, 1:]
    origin_location = node_positions[origin_index, 1:]
    # loop throught the nodes and y rotation axes
    for index, node in enumerate(node_positions):
        # vector from origin to node position
        v_op = node[1:] - origin_location
        if np.linalg.norm(v_op) > 0:  # skip the origin
            for col, axis in enumerate(axes):
                # calculate local rotation unit vector
                v_r = np.cross(axis, v_op)
                if np.linalg.norm(v_r) > 0.0:
                    v_r = v_r / np.linalg.norm(v_r)
                    d_vect = np.cross(v_r, axis)  # this should already be a unit vector
                    dist_from_axis = v_op @ d_vect / np.linalg.norm(d_vect)
                    # calculate local rotation at the centroid
                    rot = np.arctan(
                        ((all_disps[index, 1:] - origin_displacements) @ v_r)
                        / dist_from_axis
                    )
                    all_rots[col].append(rot)
                else:
                    print(
                        f"Skip point [{str(node)}] on the local rotation axis {str(axis)}."
                    )

    return all_rots


def _rotate_vector(angle, starting, axis):
    """Rotate a vector about an axis by an angle in degrees."""

    r = Rotation.from_rotvec(angle * np.array(axis), degrees=True)
    return r.apply(starting)


def _make_analysis_folder(inputs, inputs_folder, outputs_folder, name=None):

    if name is None:
        name = datetime.today().strftime("%Y-%m-%d-%H-%M-%S")  # basic timestamp
    else:
        assert isinstance(
            name, type("")
        ), "name is the output folder name and should be a string."
        # should also check that no special characters exist in name

    # create the folder and copy all inputs into it
    path = Path(outputs_folder, name)
    try:
        os.mkdir(path)
        for input_file in inputs:
            if "." not in input_file:
                input_file += ".inp"  # add extension for calculix files
            shutil.copy2(Path(inputs_folder, input_file), path)
    except OSError as error:
        print(error)

    return path


def _get_from_inp(file, keyword: str = None):
    with open(file, "r", encoding="utf-8") as f:
        contents = f.readlines()

    # extract data
    read_flag = False
    data = []
    for line in contents:
        if "*" in line and keyword in line:
            read_flag = True
            continue
        if not "**" in line and "*" in line and not keyword in line:
            break
        if read_flag:
            data.append(line.strip().rstrip(",").split(","))

    return data


def _write_nodes_file(readf, writef):
    # read the node ids from the node set file
    node_ids = np.array(
        _get_from_inp(file=readf["node_ids"], keyword="NSET,NSET="), dtype=float
    )
    # read the node locations from the mesh file in the node ids order
    all_nodes = np.array(
        _get_from_inp(file=readf["nodes"], keyword="NODE, NSET=Nall"), dtype=float
    )
    # make sure the order of the node positions is the same as the node ids,
    # as precice input forces will be read in the node id order
    sorter = np.argsort(all_nodes[:, 0])
    indices = sorter[np.searchsorted(all_nodes[:, 0], node_ids, sorter=sorter)]
    nodes = all_nodes[indices.flatten(), 1:]
    # write node ids to file
    write_to_file(file=writef["node_ids"], data=node_ids)
    # write node positions to file
    write_to_file(file=writef["nodes"], data=nodes)
    return node_ids, nodes


def _write_zero_forces_dummy(nb_nodes, writef):
    """Initialises the forces applied to the FEM to zero."""
    data = np.zeros((nb_nodes, 3))
    write_to_file(writef, data)


def _write_force_input(readf, writef):
    # read force data from precice output file
    forces = read_from_file(readf["forces"])
    node_ids = read_from_file(readf["node_ids"])
    # generate a string with the force input cards for CCX
    cards = [
        f"{int(node)},1,{force[0]:.2E}\n{int(node)},2,{force[1]:.2E}\n{int(node)},3,{force[2]:.2E}\n"
        for node, force in zip(node_ids, forces)
    ]
    # write the string to file
    with open(writef, "w", encoding="utf-8") as f:
        f.write("".join(cards))

    return forces


def _write_displacement_output(readf, writef):
    # all_disp : [nodeid, vx, vy, vz]
    all_disps = _get_from_dat(readf, data="displacements")
    write_to_file(writef, all_disps[:, 1:])


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


def _plot_forces_ditribution(nodes, forces):
    """Plot the applied force distribution in 3D."""

    zscalingf = np.max(forces[:, 2])
    ax = plt.figure().add_subplot(projection="3d")
    ax.set_xlim3d(0, np.max(nodes[:, 0]))
    ax.set_ylim3d(0, np.max(nodes[:, 1]))
    ax.set_zlim3d(0, zscalingf * 2)
    ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2], marker="o")
    ax.quiver(
        nodes[:, 0],
        nodes[:, 1],
        nodes[:, 2],
        forces[:, 0],
        forces[:, 1],
        forces[:, 2],
        normalize=False,
        arrow_length_ratio=0.0,
        color="g",
    )
    plt.show()


if __name__ == "__main__":
    main(INPUTS[2])
