"""
Parametric shell half-wing model geometry creation and meshing in CGX, followed by nastran normal modes (SOL103) analysis. 
The script assumes that CGX and nastran are installed and working. Adjust the executable calls as needed in LOCAL_EXECUTES.  
"""

# import external libraries
import warnings
import numpy as np
import csv
import subprocess
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation


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
                "elements": "EL",
                "orientation": "ORI_0",
            },
            {
                "id": "p_90",
                "thickness": 0.0002,
                "elements": "EL",
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
]

LOCAL_EXECUTES = {
    "CGX": "wsl /usr/local/bin/cgx_2.15",
    "CALCULIX": "wsl /usr/local/bin/ccx_2.19",
    "NASTRAN": "nastran",
}


def main(inputs):
    """Create a box FEM model from parametric inputs."""

    print("The parametric inputs are:")
    print(inputs)

    geometry = get_geometry(inputs, plot_flag=False)
    infile = get_CGX_input_file(geometry, inputs)
    print(f"Created cgx input file {infile}.")

    execute_CGX(infile)
    print(f"Created analysis input files with CGX.")

    if "composite_layup" in inputs:
        get_composite_properties_input(inputs)

    # run the FEM model analysis
    execute_fea(inputs["analysis_file"], inputs["fea_solver"])
    print("Executed FEM analysis.")

    # recover the analysis results
    outputs = get_fea_outputs(
        file=inputs["analysis_file"],
        solver=inputs["fea_solver"],
        mesh_file=inputs["mesh_file"],
    )
    print(outputs)

    print("End main process.\n")
    return outputs


def get_geometry(inputs, plot_flag=False):
    """Translate parameters into geometry description that CGX can understand,
    that's points, lines and surfaces."""

    aerofoil = _get_aerofoil_from_file(inputs["airfoil_csv_file"], plot_flag=plot_flag)
    points, seqa = _get_CGX_points_3D(aerofoil, inputs["chord"], inputs["span"])
    lines = _get_CGX_lines_3D(
        seqa, nele_foil=inputs["nele_foil"], nele_span=inputs["nele_span"]
    )
    surfs = _get_CGX_surfs_3D()

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

    fix_lines = [0]  # constrain the root of the wing
    loaded_lines = [1]  # apply a shear load at the tip of the wing

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
    subprocess.run(
        LOCAL_EXECUTES["CGX"] + " -bg " + infile,
        shell=True,
        check=True,
        capture_output=True,
    )


def execute_fea(file, solver):
    """ Run CGX with the batch input file to generate the mesh output files."""
    subprocess.run(
        LOCAL_EXECUTES[solver] + " " + file,
        shell=True,
        check=True,
        capture_output=True,
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


def _get_aerofoil_from_file(file, plot_flag=True):
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

    if plot_flag == 1:
        plt.plot(coordinates[:, 0], coordinates[:, 1], "-xr")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(name)
        plt.show()

    return dict(name=name, coordinates=coordinates)


def _get_CGX_points_3D(aerofoil, chord, span):
    """This function generates the CGX input file points and point sequences."""

    seqa = []

    x = aerofoil["coordinates"][:, 0] * chord
    z = aerofoil["coordinates"][:, 1] * chord

    # wing root
    y_root = np.zeros(x.size)

    # wing tip
    y_tip = np.ones(x.size) * span

    points = np.vstack([x, y_root, z]).T
    points = np.append(points, np.vstack([x, y_tip, z]).T, axis=0)

    # SEQA entries list all point IDS on the spline except for the
    # start and end points
    indices = np.arange(1, x.size - 1)
    seqa.append(indices)
    seqa.append(indices + x.size)

    return points, seqa


def _get_CGX_lines_3D(seqa, nele_foil=20, nele_span=40):
    """This function creates the aerofoil section splines and the spanwise bounding lines in CGX."""

    nele_multiplier = 2  # =2 to account for quadratic elements
    lines = []

    # aerofoil lines
    for id, seq in enumerate(seqa):
        lines.append([seq[0] - 1, seq[-1] + 1, id, nele_foil * nele_multiplier])

    # spanwise lines at trailing edge
    lines.append([seqa[0][0] - 1, seqa[1][0] - 1, nele_span * nele_multiplier])
    lines.append([seqa[0][-1] + 1, seqa[1][-1] + 1, nele_span * nele_multiplier])

    return lines


def _get_CGX_surfs_3D():

    surfaces = []
    # wing outer surface
    surfaces.append([0, 3, -1, -2])

    return surfaces


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
            line_end = " = \n" if ii + 8 <= len(points) else "\n"
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
                    if np.sign(line) > 0
                    else f"- L{np.abs(line):05d}"
                    for line in surf
                ]
            )
            + "\n"
        )

    commands.append("# =============== \n")
    # SPC and load sets
    commands.append(
        f"SETA SPC l " + " ".join([f"L{line:05d}" for line in fix_lines]) + "\n"
    )
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
    commands.append("comp SPC d\n")
    commands.append("comp LAST d\n")
    commands.append(f"send SPC {solver} spc 123456\n")
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
            f"{props['thickness']:6f},,{props['elements']},{props['orientation']}\n"
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
    main(INPUTS[0])
