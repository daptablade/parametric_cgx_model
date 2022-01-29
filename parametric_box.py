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
"""

# import external libraries
from audioop import reverse
from math import ceil
import warnings
import numpy as np
import csv
import subprocess
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation

# set these execution paths to local binaries as available
LOCAL_EXECUTES = {
    "CGX": "wsl /usr/local/bin/cgx_2.15",
    "CALCULIX": "wsl /usr/local/bin/ccx_2.19",
    "NASTRAN": "nastran",  # set to None if not installed
}

PLOT_FLAG = False

# parameters from up-stream processes and analyses
INPUTS = [
    {
        "span": 2.0,
        "chord": 0.2,
        "airfoil_csv_file": "naca0012.csv",
        "analysis_file": "normal_modes.bdf",
        "nele_foil": 20,
        "nele_span": 40,
        "node_merge_tol": 0.002,
        "cgx_ele_type": 9,  # 9: S4, 10: S8 (linear or quadratic elements)
        "cgx_solver": "nas",  # or "nas"
        "fea_solver": "NASTRAN",  # or "NASTRAN"
        "mesh_file": "all.bdf",
    },
    {
        "span": 2.0,
        "chord": 0.2,
        "airfoil_csv_file": "naca0012.csv",
        "analysis_file": "ccx_static_tip_shear",  # specify without file extension for CALCULIX
        "nele_foil": 20,
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
        "composite_layup": (
            ["p_90"] + ["p_0"] * 3 + ["p_90"] * 2 + ["p_0"] * 3 + ["p_90"]
        ),
        "composite_props_file": "composite_shell.inp",
        "mesh_file": "all.msh",
    },
    {
        "span": [0.02, 1.96, 0.02],
        "chord": [0.38, 0.38, 0.38],
        "filled_sections_flags": [True, False, True],
        "airfoil_csv_file": "naca0012.csv",
        "airfoil_cut_chord_percentages": [5, 95],
        "analysis_file": "ccx_normal_modes",  # specify without file extension for CALCULIX
        "nele_foil": [4, 10, 4, 10, 4],
        "nele_span": [4, 40, 4],
        "node_merge_tol": 0.0002,
        "cgx_ele_type": 10,  # 9: S4, 10: S8 (linear or quadratic elements)
        "cgx_solver": "abq",  # or "nas"
        "fea_solver": "CALCULIX",  # or "NASTRAN"
        "composite_plies": [
            {
                "id": "glass_0",
                "thickness": 0.0005,
                "material": "Carbon_ply",
                "orientation": "ORI_0",
            },
        ],
        "orientations": [{"id": "ORI_0", "1": [0.0, 1.0, 0.0], "2": [-1.0, 0.0, 0.0]}],
        "composite_layup": (["glass_0"]),
        "composite_props_file": "composite_shell.inp",
        "mesh_file": "all.msh",
    },
]


def main(inputs):
    """Create a box FEM model from parametric inputs."""

    print("The parametric inputs are:")
    print(inputs)

    geometry = get_geometry(inputs, plot_flag=PLOT_FLAG)
    infile = get_CGX_input_file(geometry, inputs)
    print(f"Created cgx input file {infile}.")

    execute_CGX(infile)
    print(f"Created analysis input files with CGX.")

    if "composite_layup" in inputs:
        get_composite_properties_input(inputs)

    # run the FEM model analysis
    execute_fea(inputs["analysis_file"], inputs["fea_solver"])
    print("Executed FEM analysis.")

    # # recover the analysis results
    # outputs = get_fea_outputs(
    #     file=inputs["analysis_file"],
    #     solver=inputs["fea_solver"],
    #     mesh_file=inputs["mesh_file"],
    # )
    # print(outputs)

    print("End main process.\n")
    # return outputs


def get_geometry(inputs, plot_flag=False):
    """Translate parameters into geometry description that CGX can understand,
    that's points, lines and surfaces."""

    if "airfoil_cut_chord_percentages" not in inputs:
        inputs["airfoil_cut_chord_percentages"] = None

    aerofoil = _get_aerofoil_from_file(
        inputs["airfoil_csv_file"],
        plot_flag=plot_flag,
        splitpc=inputs["airfoil_cut_chord_percentages"],
    )
    points, seqa, split_points = _get_CGX_points_3D(
        aerofoil, inputs["chord"], inputs["span"]
    )
    lines, surfs = _get_CGX_lines_3D(
        seqa,
        nele_foil=inputs["nele_foil"],
        nele_span=inputs["nele_span"],
        split_points=split_points,
    )

    return {
        "aerofoil": aerofoil,
        "points": points,
        "point_seqa": seqa,
        "lines": lines,
        "surfaces": surfs,
    }


def get_CGX_input_file(geometry, inputs):
    """Write CGX batch commands to file."""

    fdb_geom_file = "cgx_infile.fdb"

    fix_lines = None
    # [0]  # constrain the root of the wing
    loaded_lines = None
    # [1]  # apply a shear load at the tip of the wing

    # create string of all input commands
    cgx_commands = _get_commands(
        geometry,
        fix_lines,
        loaded_lines,
        merge_tol=inputs["node_merge_tol"],
        cgx_ele_type=inputs["cgx_ele_type"],
        solver=inputs["cgx_solver"],
    )

    # write string of commands to file
    with open(fdb_geom_file, "w") as f:
        f.write("".join(cgx_commands))

    return fdb_geom_file


def get_composite_properties_input(inputs):
    """write an FEA input file with the composite properties."""

    if inputs["fea_solver"] == "CALCULIX":

        # check and update the element type in the mesh input file
        str_find = "*ELEMENT, TYPE=S8,"
        str_replace = "*ELEMENT, TYPE=S8R,"
        _file_find_replace(
            file=inputs["mesh_file"], find=str_find, replace_with=str_replace
        )

        # get input file cards for this solver
        ccx_commands = _get_ccx_composite_shell_props(
            plies=inputs["composite_plies"],
            orientations=inputs["orientations"],
            layup=inputs["composite_layup"],
        )
    else:
        warnings.warn(
            f"Composite properties not implemented for solver option {inputs['fea_solver']}."
        )

    # write string of commands to file
    with open(inputs["composite_props_file"], "w") as f:
        f.write("".join(ccx_commands))


def execute_CGX(infile):
    """ Run CGX with the batch input file to generate the mesh output files."""
    if LOCAL_EXECUTES["CGX"]:
        subprocess.run(
            LOCAL_EXECUTES["CGX"] + " -bg " + infile,
            shell=True,
            check=True,
            capture_output=True,
        )
    else:
        raise ValueError("Need to specify an execution path for CalculiX GraphiX.")


def execute_fea(file, solver):
    """ Run CGX with the batch input file to generate the mesh output files."""

    if LOCAL_EXECUTES[solver]:
        subprocess.run(
            LOCAL_EXECUTES[solver] + " " + file,
            shell=True,
            check=True,
            capture_output=True,
        )
    else:
        warnings.warn(
            f"{solver} execution path has not been set. Model analysis is skipped."
        )


def get_fea_outputs(file, solver, mesh_file):
    """ Recover the analysis outputs and process them for plotting."""

    if solver == "CALCULIX":
        # read file and recover node displacements
        output_file = str(file) + ".dat"
        # all_disp : [nodeid, vx, vy, vz]
        all_disps = _get_from_dat(output_file, data="displacements")

        # average the deflections
        v_mean = np.zeros((3))
        for disp in range(3):
            v_mean[disp] = np.average(all_disps[:, disp + 1])

        # calculate the average rotations
        rotations_mean = _get_average_rotation(all_disps=all_disps, mesh_file=mesh_file)

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


def _get_aerofoil_from_file(file, plot_flag=True, splitpc=None, pt_offset=4):
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
    with open(file, mode="r") as infile:
        reader = csv.reader(infile, skipinitialspace=True)
        for row in reader:
            airfoil.append(row)
    name = airfoil[0]
    coordinates = np.array([string[0].split() for string in airfoil[1:]], dtype=float)

    # replace the last coordinate to close the airfoil at the trailing-edge
    coordinates[-1] = coordinates[0]

    splits = []
    if splitpc:
        # check that there are enough points to split the section
        min_points = 100
        if len(coordinates) < min_points:
            raise ValueError(
                f"The parameter 'airfoil_cut_chord_percentages' requires at least {min_points:d} airfoil spline points in 'airfoil_csv_file'"
            )

        # re-order the pc from TE to LE
        splitpc.sort(reverse=True)

        # we assume that there is a [0.0, 0.0] point in the airfoil
        LE_index = np.where(coordinates[:, 0] == 0.0)[0][0]
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
                        f"Values {splitpc[split_number-1]} and {split:f} in 'airfoil_cut_chord_percentages' are too close together."
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

    return dict(name=name, coordinates=coordinates, splits=splits)


def _get_CGX_points_3D(aerofoil, chord, span):
    """This function generates the CGX input file points and point sequences."""

    if not isinstance(span, list):
        span = [span]
    if not isinstance(chord, list):
        chord = [chord]

    def simple_wing(aerofoil, chord, span):
        seqa = []
        starting_y = 0
        pt_counter = 0
        for section_index, length_y in enumerate(span):

            x = aerofoil["coordinates"][:, 0] * chord[section_index]
            z = aerofoil["coordinates"][:, 1] * chord[section_index]

            if section_index == 0:  # only needed at the root of the wing
                y_root = np.ones(x.size) * starting_y
                points = np.vstack([x, y_root, z]).T

            # section tip
            y_tip = np.ones(x.size) * (starting_y + length_y)
            points = np.append(points, np.vstack([x, y_tip, z]).T, axis=0)

            # SEQA entries list all point IDS on the spline except for the
            # start and end points
            indices = np.arange(pt_counter + 1, pt_counter + x.size - 1)
            seqa.append(indices)
            if section_index == 0:  # only needed at the root of the wing
                seqa.append(indices + x.size)

            starting_y += length_y
            pt_counter = seqa[-1][-1] + 2

            return points, seqa

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

            def airfoil_SEQA(pt_counter, seqa, all_split_points):
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
                        pt_counter + aerofoil["splits"][-1]["bot"],
                    )
                )
                # SEQA for first airfoil bot spline
                pt = x.size - 1
                bot_seqa = []
                for split in aerofoil["splits"]:
                    indices = np.flipud(
                        np.arange(pt_counter + pt - 1, pt_counter + split["bot"], -1)
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
                return seqa, all_split_points

            seqa, split_points = airfoil_SEQA(
                pt_counter=pt_counter, seqa=seqa, all_split_points=split_points
            )
            if section_index == 0:  # only needed at the root of the wing
                seqa, split_points = airfoil_SEQA(
                    pt_counter=pt_counter + x.size,
                    seqa=seqa,
                    all_split_points=split_points,
                )

            starting_y += length_y
            pt_counter = seqa[-1][-1] + 2

        return points, seqa, split_points

    if "splits" not in aerofoil:
        points, seqa = simple_wing(aerofoil, chord, span)
        split_points = None
    else:
        points, seqa, split_points = wing_with_splits(aerofoil, chord, span)

    return points, seqa, split_points


def _get_CGX_lines_3D(
    seqa, nele_foil=20, nele_span=40, nele_split=4, split_points=None
):
    """This function creates the aerofoil section splines and the spanwise bounding lines in CGX."""

    nele_multiplier = 2  # =2 to account for quadratic elements
    lines = []
    aero_surfaces = []
    rib_surfaces = []

    if not isinstance(nele_span, list):
        nele_span = [nele_span]

    if not isinstance(nele_foil, list):
        nele_foil = [nele_foil]

    splits = 0
    seqas_per_aerofoil = 1
    aerofoils = len(seqa)
    if isinstance(split_points, np.ndarray):
        splits = split_points.shape[0]
        seqas_per_aerofoil = splits * 2 + 1
        aerofoils = int(len(seqa) / seqas_per_aerofoil)

    airfoil_index = 0
    lcounter = 0
    for id, seq in enumerate(seqa):

        # aerofoil lines
        lines.append(
            [
                seq[0] - 1,
                seq[-1] + 1,
                id,
                nele_foil[id % seqas_per_aerofoil] * nele_multiplier,
            ]
        )
        lcounter += 1

        if (id + 1) % seqas_per_aerofoil == 0:

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
            if (id + 1) / seqas_per_aerofoil < aerofoils:

                for te_line_inc in range(seqas_per_aerofoil + 1):
                    if te_line_inc < seqas_per_aerofoil:
                        start_id = id + 1 - seqas_per_aerofoil + te_line_inc
                        end_id = id + 1 + te_line_inc
                        side = 0
                        pt_offset = -1
                    else:
                        start_id = id - seqas_per_aerofoil + te_line_inc
                        end_id = id + te_line_inc
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

    # surfaces
    surfaces = aero_surfaces + rib_surfaces

    return lines, surfaces


def _get_commands(
    geometry, fix_lines, loaded_lines, merge_tol=0.001, cgx_ele_type=10, solver="abq"
):

    commands = []

    # points
    for id, point in enumerate(geometry["points"]):
        commands.append(f"PNT P{id:05d} {point[0]:e} {point[1]:e} {point[2]:e}\n")

    commands.append("# =============== \n")
    # point sequences
    for id, points in enumerate(geometry["point_seqa"]):
        commands.append(f"SEQA A{id:05d} pnt ")
        for ii in range(0, len(points), 8):
            line_end = " = \n" if ii + 8 < len(points) else "\n"
            commands.append(
                " ".join([f"P{point:05d}" for point in points[ii : ii + 8]]) + line_end
            )

    commands.append("# =============== \n")
    # lines
    for id, line in enumerate(geometry["lines"]):
        if len(line) == 3:  # straight line
            commands.append(
                f"LINE L{id:05d} P{line[0]:05d} P{line[1]:05d} {line[2]:d} \n"
            )
        elif len(line) == 4:  # spline
            commands.append(
                f"LINE L{id:05d} P{line[0]:05d} P{line[1]:05d} A{line[2]:05d} {line[3]:d} \n"
            )

    commands.append("# =============== \n")
    # surfaces
    for id, surf in enumerate(geometry["surfaces"]):
        commands.append(
            f"GSUR V{id:05d} + BLEND "
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

    commands.append("# =============== \n")
    # SPC and load sets
    if fix_lines:
        commands.append(
            f"SETA SPC l " + " ".join([f"L{line:05d}" for line in fix_lines]) + "\n"
        )
    if loaded_lines:
        commands.append(
            f"SETA LAST l " + " ".join([f"L{line:05d}" for line in loaded_lines]) + "\n"
        )

    commands.append("# =============== \n")
    # surface meshes
    for id, _ in enumerate(geometry["surfaces"]):
        commands.append(f"MSHP V{id:05d} s {cgx_ele_type:d} 0 1.000000e+00\n")

    commands.append("# =============== \n")
    # custom export statement
    commands.append("mesh all\n")
    commands.append(f"merg n all {merge_tol:6f} 'nolock'\n")
    commands.append("comp nodes d\n")
    if fix_lines:
        commands.append("comp SPC d\n")
        commands.append(f"send SPC {solver} spc 123456\n")
    if loaded_lines:
        commands.append("comp LAST d\n")
        commands.append(f"send LAST {solver} names\n")
    commands.append(f"send all {solver} \n")
    commands.append("quit\n")

    return commands


def _get_ccx_composite_shell_props(plies=None, orientations=None, layup=None):

    commands = []

    # orientation cards
    for ori in orientations:
        commands.append(f"*ORIENTATION,NAME={ori['id']}\n")
        commands.append(", ".join(str(x) for x in [*ori["1"], *ori["2"]]) + "\n")

    commands.append("** =============== \n")
    # shell property
    commands.append("*SHELL SECTION,ELSET=Eall,COMPOSITE\n")
    for ply in layup:
        props = [p for p in plies if p["id"] == ply][0]
        commands.append(
            f"{props['thickness']:6f},,{props['material']},{props['orientation']}\n"
        )

    return commands


def _file_find_replace(file, find: str, replace_with: str):
    with open(file, "r") as f:
        contents = f.readlines()

    for index, line in enumerate(contents):
        if find in line:
            contents[index] = line.replace(find, replace_with)
            print(f"Find & Replace edited file '{file}' at line {index:d}.")
            break

    with open(file, "w") as f:
        f.write("".join(contents))


def _get_from_dat(file, data: str = None):
    with open(file, "r") as f:
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

    origin_position = np.zeros((3))
    # calculate the approximate rotation from each node
    rotations = _get_rot_from_disp(all_disps, node_positions, axes, origin_position)

    # average the rotation values
    rotations_mean = np.zeros((3))
    for rot in range(3):
        rotations_mean[rot] = np.average(rotations[:, rot])

    return rotations_mean


def _get_nodes_from_inp(file, ref_nodes=None):
    with open(file, "r") as f:
        contents = f.readlines()

    # find node definitions in the file
    start = None
    end = None
    nodes = []
    for index, line in enumerate(contents):
        if "*NODE" in line:
            start = index + 1
        elif start and "*" in line:  # marks the next keyword
            end = index - 1
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


def _get_rot_from_disp(all_disps, node_positions, axes, origin_position):
    """Return an array of rotations in radians. Using small deflection assumptions."""

    all_rots = np.zeros((node_positions.shape[0], 3))
    # loop throught the nodes and y rotation axes
    for index, node in enumerate(node_positions):
        # vector from origin to node position
        v_op = node[1:] - origin_position
        for col, axis in enumerate(axes):
            # calculate local rotation unit vector
            v_r = np.cross(axis, v_op)
            v_r = v_r / np.linalg.norm(v_r)
            # calculate local rotation at the centroid
            rot = (all_disps[index, 1:] @ v_r) / np.linalg.norm(v_op)
            if np.abs(rot) >= 1:
                raise ValueError("Rotation error.")
            all_rots[index, col] = rot

    all_rots = np.arcsin(all_rots)

    return all_rots


def _rotate_vector(angle, starting, axis):
    """Rotate a vector about an axis by an angle in degrees."""

    r = Rotation.from_rotvec(angle * np.array(axis), degrees=True)
    return r.apply(starting)


if __name__ == "__main__":
    main(INPUTS[2])
