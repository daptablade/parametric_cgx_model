"""Plot precice and solver iteration output data."""

from pathlib import Path


def write_solver_output_to_file(solver, output, folder):
    """write solver iteration outputs to file."""
    with open(folder / (solver + "_iter_outputs.txt"), "w", encoding="utf-8") as f:
        for i, out in enumerate(output):
            f.write(str(i) + ": " + out + "\n")


def read_solver_outputs(solver, folder):
    """read solver iteration outputs from file."""
    with open(folder / (solver + "_iter_outputs.txt"), "r", encoding="utf-8") as f:
        contents = f.readlines()
    return contents


def post_process_solver_iters(solver, folder):
    """select the post-processing function depending on the solver."""
    output = read_solver_outputs(solver, folder)
    if solver == "SolverOne":
        _post_process_strip_theory_aero(output)
    elif solver == "SolverTwo":
        _post_process_param_box(output)
    else:
        raise ValueError(f"{solver} is not a valid post-processing option.")


def _post_process_param_box(output):
    print(output)


def _post_process_strip_theory_aero(output):
    print(output)


if __name__ == "__main__":
    post_process_solver_iters(
        solver="SolverTwo", folder=Path(__file__).parent / "outputs"
    )
