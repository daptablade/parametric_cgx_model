"""
Parametric shell half-wing model geometry creation and meshing in CGX, followed by nastran normal modes (SOL103) analysis. 
The script assumes that CGX and nastran are installed and working. Adjust the executable calls as needed CGX_EXECUTE and NASTRAN_EXECUTE.  
"""

# import external libraries
import numpy as np
import csv
import subprocess
import matplotlib.pyplot as plt


# parameters from up-stream processes and analyses
INPUTS = {
    "span": 2.0,
    "chord": 0.2,
    "airfoil_csv_file": "naca0012.csv",
    "analysis_file": "normal_modes.bdf",
    "nele_foil": 20,
    "nele_span": 40,
    "node_merge_tol": 0.002,
}

CGX_EXECUTE = "wsl /usr/local/bin/cgx_2.15"
NASTRAN_EXECUTE = "nastran"


def main(inputs):
    """Create a box FEM model from parametric inputs."""

    print("The parametric inputs are: \n")
    print(inputs)

    geometry = get_geometry(inputs)
    infile = get_CGX_input_file(geometry, inputs)
    print(f"Created cgx input file {infile}.")

    execute_CGX(infile)
    print(f"Created analysis input files with CGX.")

    # run the FEM model analysis
    execute_nastran(inputs["analysis_file"])
    print(f"Executed FEM analysis.\nEnd main process.")


def get_geometry(inputs):
    """Translate parameters into geometry description that CGX can understand,
    that's points, lines and surfaces."""

    aerofoil = _get_aerofoil_from_file(inputs["airfoil_csv_file"])
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

    # create string of all input commands
    cgx_commands = _get_commands(
        geometry, fix_lines, merge_tol=inputs["node_merge_tol"]
    )

    # write string of commands to file
    with open(fdb_geom_file, "w") as f:
        f.write("".join(cgx_commands))

    return fdb_geom_file


def execute_CGX(infile):
    """ Run CGX with the batch input file to generate the mesh output files."""
    subprocess.run(
        CGX_EXECUTE + " -bg " + infile,
        shell=True,
        check=True,
        capture_output=True,
    )


def execute_nastran(file):
    """ Run CGX with the batch input file to generate the mesh output files."""
    subprocess.run(
        NASTRAN_EXECUTE + " " + file,
        shell=True,
        check=True,
        capture_output=True,
    )


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
    coordinates = np.array(
        [string[0].split() for string in airfoil[1:]], dtype=np.float
    )

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


def _get_commands(geometry, fix_lines, merge_tol=0.001):

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
    # SPC sets
    commands.append(
        f"SETA SPC l " + " ".join([f"L{line:05d}" for line in fix_lines]) + "\n"
    )

    commands.append("# =============== \n")
    # surface meshes
    for id, _ in enumerate(geometry["surfaces"]):
        commands.append(f"MSHP V{id:05d} s 9 0 1.000000e+00\n")

    commands.append("# =============== \n")
    # custom export statement
    commands.append("mesh all\n")
    commands.append(f"merg n all {merge_tol:6f} 'nolock'\n")
    commands.append("comp nodes d\n")
    commands.append("comp SPC d\n")
    commands.append("send SPC nas spc 123456\n")
    commands.append("send all nas\n")
    commands.append("quit\n")

    return commands


if __name__ == "__main__":
    main(INPUTS)
