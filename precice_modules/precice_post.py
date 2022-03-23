"""Plot precice and solver iteration output data."""

from pathlib import Path
import ast
import local_context
from parametric_studies import _plot_study_results


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
        return _post_process_strip_theory_aero(output, folder)
    elif solver == "SolverTwo":
        return _post_process_param_box(output, folder)
    else:
        raise ValueError(f"{solver} is not a valid post-processing option.")


def _post_process_param_box(output, folder):
    parsed_output = []
    for line in output:
        if "{'Ux':" in line:
            parsed_output.append(ast.literal_eval(line))
    _plot_study_results(
        parsed_output,
        x=range(len(parsed_output)),
        y=["Uz", "Ry"],
        xlabel="iteration",
        ylabel="displacements (m) or rotation (rad)",
        saveas=str(folder / "wing_tip_displacements_convergence"),
    )
    print("Param_box converged output: ", parsed_output[-1])
    return parsed_output[-1]


def _post_process_strip_theory_aero(output, folder):
    parsed_output = []
    iterator = 0
    read_aoa = False
    read_lift = False
    for line in output:
        if "The parametric inputs are" in line:
            if iterator != int(line[0]):
                parsed_output.append(
                    {
                        "AoA": [float(a) for a in aoa],
                        "Lift": [float(a) for a in lift],
                    }
                )
                iterator = int(line[0])
            aoa = []
            lift = []
        if "Local strip angle of incidence in degrees:" in line:
            read_aoa = True
            continue
        if "Local strip lift force at the aerodynamic centre:" in line:
            read_lift = True
            continue
        if read_aoa:
            if "]" in line:
                read_aoa = False
            aoa += line.replace("[", "").replace("]", "").split()
        if read_lift:
            if "]" in line:
                read_lift = False
            lift += line.replace("[", "").replace("]", "").split()

    _plot_study_results(
        [parsed_output[-1]],
        x=range(len(parsed_output[-1]["AoA"])),
        y=["AoA"],
        xlabel="spanwise strip (root to tip)",
        ylabel="Aerofoil incidence (deg)",
        saveas=str(folder / "aero_incidence_distribution"),
    )
    _plot_study_results(
        parsed_output,
        x=range(len(parsed_output)),
        y=["AoA"],
        xlabel="iteration",
        ylabel="Aerofoil incidence (deg)",
        saveas=str(folder / "aero_incidence_convergence"),
    )
    _plot_study_results(
        [parsed_output[-1]],
        x=range(len(parsed_output[-1]["Lift"])),
        y=["Lift"],
        xlabel="spanwise strip (root to tip)",
        ylabel="Lift (N)",
        saveas=str(folder / "aero_lift_distribution"),
    )
    _plot_study_results(
        parsed_output,
        x=range(len(parsed_output)),
        y=["Lift"],
        xlabel="iteration",
        ylabel="Lift (N)",
        saveas=str(folder / "aero_lift_convergence"),
    )
    print("aero converged output: ", parsed_output[-1])
    return parsed_output[-1]


if __name__ == "__main__":
    converged_aero_out = post_process_solver_iters(
        solver="SolverOne", folder=Path(__file__).parent / "outputs"
    )
    converged_strucutures_out = post_process_solver_iters(
        solver="SolverTwo", folder=Path(__file__).parent / "outputs"
    )
