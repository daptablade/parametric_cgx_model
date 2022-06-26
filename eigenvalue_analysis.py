from scipy.sparse import csr_matrix, csc_matrix, bmat, eye
from scipy.sparse.linalg import eigsh, eigs, inv
from scipy.interpolate import interp2d
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import csv

from utils import timeit
from strip_theory_aerodynamics import AeroModel, AERO_INPUTS


PLATE_AERO = {
    "planform": {
        "A": np.array([0.0, 0.0, 0.0], dtype=float),  # LE fixed root
        "B": np.array([0.0, 0.3048, 0.0], dtype=float),  # LE tip
        "CA": 0.0762,  # chord at root
        "CB": 0.0762,  # chord at tip
    },
    "strips": 10,  # number of aero strips along the span >=1
    "root_alpha": 10,  # AoA at the wing root in degrees
    "rho": 1.225,  # air density in Pa
    "V": 5.0,  # air velocity in m/s
    "CL_alpha": 2 * np.pi,  # ideal lift curve slope
}


def read_matrix_csv(file, delimiter=" ", skip_lines_with=None, read_only=None):
    """Reads csv and strips empty entries."""

    if read_only is None:
        read_flag = True
    else:
        read_flag = False

    with open(file, newline="") as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter)
        if skip_lines_with is None and read_flag:
            data = [[float(val) for val in row if bool(val.strip())] for row in reader]
        elif skip_lines_with and read_flag:
            data = [
                [float(val) for val in row if bool(val.strip())]
                for row in reader
                if not skip_lines_with in row[0]
            ]
        elif not read_flag:
            data = []
            for row in reader:
                if not read_flag and read_only in row[0]:
                    read_flag = True
                    continue
                elif read_flag and skip_lines_with in row[0]:
                    read_flag = False
                    break
                elif read_flag:
                    data.append([float(val) for val in row if bool(val.strip())])

    return np.array(data)


@timeit
def get_matrix(data, filter_out_rows=None):

    mask = np.ones(len(data), dtype=bool)
    if filter_out_rows:
        # set the mask to false where the index is a filtered row
        frows = []
        for row in filter_out_rows:
            # get the filtered row / col indices
            index = np.sort(
                np.hstack(
                    [np.where(data[:, 0] == row)[0], np.where(data[:, 1] == row)[0]]
                )
            )
            # index = np.sort(np.where(data[:,1] == row)[0])
            if index.size > 0:
                frows += index.tolist()
        mask[frows] = False

    # data array contains (row, col, val) in each row
    m, n = np.max(data[:, :2], axis=0)
    matrix = np.zeros((int(m), int(n)), dtype=np.float64)
    for row in data[mask, :]:
        matrix[int(row[0] - 1), int(row[1] - 1)] = row[2]
    submatrix = np.where(np.any(matrix < 0, axis=1) | np.any(matrix > 0, axis=1))[0]
    mat = matrix[np.reshape(submatrix, (-1, 1)), submatrix]
    return csr_matrix(mat + np.tril(mat.T, k=-1))


def get_rows_from_dofs(dofs, lookup):
    rows = []
    for dof in dofs:
        index = np.where((lookup == dof).all(axis=1))[0] + 1
        if index.size > 0:
            rows.append(int(index))
    return rows


@timeit
def get_normal_modes(problem):
    # solve an eigenvalue problem to obtain 10 lowest eigenvalues
    # evals_small, evecs_small = eigs(
    #     problem["sti"], k=10, M=problem["mas"], sigma=0.0, which="LR"
    # )
    # # print frequencies in Hz
    # omega = np.sqrt(evals_small) / (2 * np.pi)
    # print("Eigenvalues (Hz): ", omega)

    evals_small, evecs_small = eigsh(
        problem["sti"], k=10, M=problem["mas"], sigma=0.0, which="LM"
    )
    # print frequencies in Hz
    omega = np.sqrt(evals_small) / (2 * np.pi)
    print("Eigenvalues (Hz): ", omega)

    return omega, evals_small, evecs_small


def filter_evec(node_ids, dirs, lookup, evecs):

    mask = np.ones(len(lookup), dtype=bool)
    for index, row in enumerate(lookup):
        if not row[0] in node_ids or not row[1] in dirs:
            mask[index] = False

    filtered_evecs = evecs[mask, :]
    filtered_dofs = lookup[mask, :]

    return filtered_evecs, filtered_dofs


def get_unique_evecs(filtered_evecs, filtered_dofs):

    reduced_dofs = np.unique(filtered_dofs, axis=0)
    reduced_evecs = np.zeros([len(reduced_dofs), filtered_evecs.shape[1]])
    for index, dof in enumerate(reduced_dofs):
        # get indices of repeated dofs
        indices = np.where((filtered_dofs == dof).all(axis=1))[0]
        reduced_evecs[index, :] = np.average(filtered_evecs[indices, :], axis=0)

    return reduced_evecs, reduced_dofs


def apply_boundary_conditions(reduced_evecs, reduced_dofs, fixed):

    mask = np.zeros(len(reduced_dofs), dtype=bool)
    frows = []
    for row in fixed:
        index = np.where((reduced_dofs == row).all(axis=1))[0]
        # index = np.sort(np.where(data[:,1] == row)[0])
        if index.size > 0:
            frows += index.tolist()
    mask[frows] = True

    reduced_evecs[mask, :] = 0.0

    return reduced_evecs


def plot_modes(data_sets, modes="all"):

    if modes == "all":
        # one figure per eigenvector
        modes = data_sets[0][2].shape[1]

    for mode in range(modes):
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        for data in data_sets:
            ax.scatter(data[0], data[1], data[2][:, mode], label=data[3])
            if data[3] == "AERO":
                LE_indices = int(len(data[0]) / 2)
                ax.plot3D(
                    data[0][:LE_indices],
                    data[1][:LE_indices],
                    data[2][:LE_indices, mode],
                    "grey",
                )
                ax.plot3D(
                    data[0][LE_indices:],
                    data[1][LE_indices:],
                    data[2][LE_indices:, mode],
                    "grey",
                )
        ax.set_xlabel("X, m")
        ax.set_ylabel("Y, m")
        ax.set_zlabel("Z, normalised")
        ax.legend()

    plt.show()


def get_aero_evects(
    evecs, lookup, nodes_xyz, aeronodes, fixed=None, dirs=[3], plot_flag=False
):

    # filter evecs to nodes xyz and dims only
    filtered_evecs, filtered_dofs = filter_evec(
        node_ids=nodes_xyz[:, 0].astype(int),
        dirs=dirs,
        lookup=lookup.astype(int),
        evecs=evecs,
    )

    # average duplicate dofs
    reduced_evecs, reduced_dofs = get_unique_evecs(
        filtered_evecs, filtered_dofs.astype(int)
    )

    # set boundary conditions to zero
    if not fixed is None:
        reduced_evecs = apply_boundary_conditions(reduced_evecs, reduced_dofs, fixed)

    # get xyz for all dof nodes
    X = np.zeros([len(reduced_dofs)])
    Y = np.zeros([len(reduced_dofs)])
    for index, dof in enumerate(reduced_dofs):
        node_index = np.where(nodes_xyz[:, 0] == dof[0])[0]
        X[index] = nodes_xyz[node_index, 1]
        Y[index] = nodes_xyz[node_index, 2]

    # interpolate aero eigenvectors
    aero_evecs = np.zeros([len(aeronodes), reduced_evecs.shape[1]])
    for index, evect in enumerate(reduced_evecs.transpose()):
        f = interp2d(X, Y, evect, kind="linear")
        for nodeid, node in enumerate(aeronodes):
            aero_evecs[nodeid, index] = f(node[0], node[1])

    # plot selected modes for visual check
    if plot_flag:
        plot_modes(
            [
                (X, Y, reduced_evecs, "FEM"),
                (aeronodes[:, 0], aeronodes[:, 1], aero_evecs, "AERO"),
            ]
        )

    return aero_evecs


def get_chords_from_aero_dofs(aeromodel):
    c = np.zeros([len(aeromodel.le_nodes), 1])
    for ii in range(len(aeromodel.le_nodes)):
        c[ii] = aeromodel.te_nodes[ii, 0] - aeromodel.le_nodes[ii, 0]
    return c


def get_lift_factor_from_aero_dofs(aeromodel, span):
    aw = np.zeros([len(aeromodel.le_nodes), 1])
    for ii in range(len(aeromodel.le_nodes)):
        aw[ii] = 2 * np.pi * (1 - (aeromodel.le_nodes[ii, 1] / span) ** 3)
    return aw


def get_phi(aero_evecs, option):

    nnodes = int(aero_evecs.shape[0] / 2)
    phi = np.zeros([nnodes, aero_evecs.shape[1]])
    for index, col in enumerate(aero_evecs.T):
        if option == "LEplusTE":
            phi[:, index] = col[:nnodes] + col[nnodes:]
        elif option == "LEminusTE":
            phi[:, index] = col[:nnodes] - col[nnodes:]
        else:
            raise ValueError(f"Option {option} is not recognised.")
    return phi


def get_aero_component_matrices(aeromodel, aero_evecs):

    # get component matrices
    span = aeromodel._point_b[1]
    nstrips = aeromodel.aero_inputs["strips"]
    c = get_chords_from_aero_dofs(aeromodel)  # vector of length dof/2
    aw = get_lift_factor_from_aero_dofs(aeromodel, span)  # vector of length dof/2
    dy_strips = (
        np.vstack([np.array([0.5]), np.ones([nstrips - 1, 1]), np.array([0.5])])
        * span
        / (nstrips + 1)
    )  # vector of length dof/2

    # matrices of dof/2 * n eigenvectors
    phi_LEplusTE = get_phi(aero_evecs, option="LEplusTE")
    phi_LEminusTE = get_phi(aero_evecs, option="LEminusTE")

    return c, aw, dy_strips, phi_LEplusTE, phi_LEminusTE


def get_aero_damping(aeromodel, aero_evecs, Mxi_dot=-1.2, e=0.25):
    c, aw, dy_strips, phi_LEplusTE, phi_LEminusTE = get_aero_component_matrices(
        aeromodel, aero_evecs
    )
    mat_size = aero_evecs.shape[1]
    B_mat = np.zeros([mat_size, mat_size], dtype=np.float64)
    for index in range(mat_size):
        B_mat[index, :] = (
            (1 / 8)
            * phi_LEplusTE[:, index].transpose()
            @ np.diag((c * aw * dy_strips)[:, 0])
            @ phi_LEplusTE
            - (e / 4)
            * phi_LEminusTE[:, index].transpose()
            @ np.diag((c * aw * dy_strips)[:, 0])
            @ phi_LEplusTE
            - (Mxi_dot / 8)
            * phi_LEminusTE[:, index].transpose()
            @ np.diag((c * dy_strips)[:, 0])
            @ phi_LEminusTE
        )

    return csr_matrix(B_mat)


def get_aero_stiffness(aeromodel, aero_evecs, e=0.25):
    _, aw, dy_strips, phi_LEplusTE, phi_LEminusTE = get_aero_component_matrices(
        aeromodel, aero_evecs
    )
    mat_size = aero_evecs.shape[1]
    C_mat = np.zeros([mat_size, mat_size], dtype=np.float64)
    for index in range(mat_size):
        C_mat[index, :] = (1 / 4) * phi_LEplusTE[:, index].transpose() @ np.diag(
            (aw * dy_strips)[:, 0]
        ) @ phi_LEminusTE - (e / 2) * phi_LEminusTE[:, index].transpose() @ np.diag(
            (aw * dy_strips)[:, 0]
        ) @ phi_LEminusTE

    return csr_matrix(C_mat)


def get_system_matrix(problem, rho, V):

    M = problem["M"]
    K = problem["K"]
    B = problem["B"]
    C = problem["C"]
    Q = bmat(
        [
            [None, eye(M.shape[0], dtype=np.float64)],
            [-inv(M) @ (K + rho * V**2 * C), -inv(M) @ (rho * V * B)],
        ]
    )

    return Q


@timeit
def get_complex_aero_modes(problem, rho, V):

    Q = get_system_matrix(problem, rho, V)

    evals_small, evecs_small = eigs(Q, k=10, sigma=0.0, which="LM")
    # damping ratios
    # print frequencies in Hz
    V_omega = np.abs(evals_small) / (2 * np.pi)
    print("Eigenvalues (Hz): ", V_omega)

    return V_omega


def main(file, folder, aero_inputs=None, box_inputs=None):

    # dofs that are constrained
    csv_file = Path(folder, "SPC_123456.bou")
    fixed_dofs = read_matrix_csv(csv_file, delimiter=",", skip_lines_with="**")

    # load row_dof lookup matrix, where each row matches K, M matrix rows and corresponds
    # to node.dir combination
    csv_file = Path(folder, file + "." + "dof")
    lookup = read_matrix_csv(csv_file, delimiter=".")

    # read shell node locations from input file
    csv_file = Path(folder, "all.msh")
    nodes_xyz = read_matrix_csv(
        csv_file, delimiter=",", skip_lines_with="*", read_only="*NODE"
    )
    # filter nodes for aero shape interpolation (e.g. only upper surface on box)
    nodes_xyz = nodes_xyz[nodes_xyz[:, 3] >= 0.0]

    problem = {}
    for matrix in ["mas", "sti"]:
        csv_file = Path(folder, file + "." + matrix)
        data = read_matrix_csv(csv_file)
        problem[matrix] = get_matrix(data)  # , filter_out_rows=fixed_dofs_rows)

    # check that the stiffness matrix is positive definite
    if not all(eigsh(problem["sti"])[0] > 0):
        print("Warning: Stiffness matrix is not positive definite.")
    if not all(eigsh(problem["mas"])[0] > 0):
        print("Warning: Mass matrix is not positive definite.")

    # normal modes analysis
    omega, evals, evecs = get_normal_modes(problem)

    if aero_inputs:
        # get modal mass and stiffness matrices
        problem["M"] = csc_matrix(evecs.T @ problem["mas"] @ evecs)
        problem["K"] = csc_matrix(evecs.T @ problem["sti"] @ evecs)

        # instantiate aero model
        aeromodel = AeroModel(aero_inputs, box_inputs)

        # define the corner nodes of the aerodynamic planform
        aeromodel.get_planform()

        # calculate the internal aero node positions at the leading and trailing edges
        aeromodel.get_all_nodes()
        aeronodes = np.vstack([aeromodel.le_nodes, aeromodel.te_nodes])

        # filter normal mode eigenvector DOFs and reduce from 3D element nodes to shell nodes
        aero_evecs = get_aero_evects(
            evecs, lookup, nodes_xyz, aeronodes, fixed=fixed_dofs, plot_flag=False
        )

        # get aero damping and stiffness matrices B and C
        problem["B"] = get_aero_damping(aeromodel, aero_evecs, Mxi_dot=-1.2, e=0.25)
        problem["C"] = get_aero_stiffness(aeromodel, aero_evecs, e=0.25)

        # flutter / divergence analysis
        V_omega = get_complex_aero_modes(
            problem, rho=aero_inputs["rho"], V=aero_inputs["V"]
        )

    return omega

    # results = eig(problem["sti"], b=problem["mas"])
    # print(np.sqrt(np.sort(results[0])) / (2 * np.pi))


if __name__ == "__main__":

    # SIMPLE PLATE normal modes checks
    freq_scipy = main(
        file="ccx_normal_modes_matout", folder="outputs/test_mat_out/test_simple_plate"
    )
    freq_ccx = np.array(
        [
            0.1112296e02,
            0.3969891e02,
            0.7415313e02,
            0.1408884e03,
            0.3089473e03,
            0.4026458e03,
            0.7748749e03,
            0.2052274e04,
            0.2144465e04,
            0.3231235e04,
        ]
    )
    assert np.allclose(freq_scipy, freq_ccx, rtol=1e-3)

    # SIMPLE PLATE flutter analysis check
    freq_scipy = main(
        file="ccx_normal_modes_matout",
        folder="outputs/test_mat_out/test_simple_plate",
        aero_inputs=PLATE_AERO,
    )

    # box model normal modes checks
    freq_scipy = main(file="ccx_normal_modes_matout", folder="outputs/test_mat_out")
    freq_ccx = np.array(
        [
            0.8692880e01,
            0.4375812e02,
            0.5737775e02,
            0.6076665e02,
            0.9949578e02,
            0.1538641e03,
            0.1741874e03,
            0.2012603e03,
            0.2374603e03,
            0.2688630e03,
        ]
    )
    assert np.allclose(freq_scipy, freq_ccx, rtol=1e-3)
