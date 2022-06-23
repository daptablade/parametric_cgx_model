from scipy.sparse import coo_matrix, bsr_matrix, dok_matrix
from scipy.sparse.linalg import eigsh, eigs
from scipy.linalg import eig, eigh
import numpy as np

from pathlib import Path
import csv


def read_matrix_csv(file, delimiter=" ", skip_lines_with=None):
    """Reads csv and strips empty entries."""

    with open(file, newline="") as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter)
        if skip_lines_with is None:
            data = [[float(val) for val in row if bool(val.strip())] for row in reader]
        else:
            data = [
                [float(val) for val in row if bool(val.strip())]
                for row in reader
                if not skip_lines_with in row[0]
            ]
    return np.array(data)


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
    return mat + np.tril(mat.T, k=-1)


# def node_dir_2_index(node_dir, dof_per_node=3):
#     return int((node_dir[0] - 1) * dof_per_node + node_dir[1] - 1)


def get_rows_from_dofs(dofs, lookup):
    rows = []
    for dof in dofs:
        index = np.where((lookup == dof).all(axis=1))[0] + 1
        if index.size > 0:
            rows.append(int(index))
    return rows


def main(
    file="ccx_normal_modes_matout", folder="outputs/test_mat_out/test_simple_plate"
):

    # dofs that are constrained
    csv_file = Path(folder, "SPC_123456.bou")
    fixed_dofs = read_matrix_csv(csv_file, delimiter=",", skip_lines_with="**")

    # row entries for dofs that are constrained
    csv_file = Path(folder, file + "." + "dof")
    lookup = read_matrix_csv(csv_file, delimiter=".")
    fixed_dofs_rows = get_rows_from_dofs(
        dofs=fixed_dofs.astype(int), lookup=lookup.astype(int)
    )

    problem = {}
    for matrix in ["mas", "sti"]:
        csv_file = Path(folder, file + "." + matrix)
        data = read_matrix_csv(csv_file)
        problem[matrix] = get_matrix(data)  # , filter_out_rows=fixed_dofs_rows)

    # check that the stiffness matrix is positive definite
    if not all(eig(problem["sti"])[0] > 0):
        print("Warning: Stiffness matrix is not positive definite.")
    if not all(eig(problem["mas"])[0] > 0):
        print("Warning: Mass matrix is not positive definite.")

    # solve an eigenvalue problem to obtain 10 lowest eigenvalues
    evals_small, evecs_small = eigs(
        problem["sti"], k=10, M=problem["mas"], sigma=0.0, which="LR"
    )
    # print frequencies in Hz
    omega = np.sqrt(evals_small) / (2 * np.pi)
    print("Eigenvalues (Hz): ", omega)

    evals_small, evecs_small = eigsh(
        problem["sti"], k=10, M=problem["mas"], sigma=0.0, which="LM"
    )
    # print frequencies in Hz
    omega = np.sqrt(evals_small) / (2 * np.pi)
    print("Eigenvalues (Hz): ", omega)

    # results = eig(problem["sti"], b=problem["mas"])
    # print(np.sqrt(np.sort(results[0])) / (2 * np.pi))


if __name__ == "__main__":
    main()
