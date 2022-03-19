"""Execute the python script 'parametric_studies.py' and inspect outputs:
this script executes multiple iterations of the main() function from parametric_box.py
with incremental input changes. Note: Calculix CrunchiX is required for this script and
the exectuion path needs to be set in LOCAL_EXECUTES in parametric_studies.py.

Choose between: 
>> parametric_study_material(INPUTS[1]) for a composite material design study
>> parametric_study_mesh(INPUTS[1]) for a mesh refinement study
"""

from shutil import copy2
import copy
import numpy as np
import matplotlib.pyplot as plt

from parametric_box import INPUTS, main, _rotate_vector, _file_find_replace


def parametric_study_material(inputs):
    """Parametric study exploring the effect of material properties on model deflections."""

    study_results = []
    params = [-30, -10, 0, 10, 30]  # rotation_angles, deg

    for iteration, param_value in enumerate(params):
        print(f"Parametric study iteration {iteration}.")
        inputs_copy = copy.deepcopy(inputs)

        inputs_copy["orientations"][0]["1"] = _rotate_vector(
            angle=param_value,
            starting=inputs_copy["orientations"][0]["1"],
            axis=[0.0, 0.0, 1.0],
        )

        output = main(inputs_copy)
        output["inputs"] = inputs_copy
        study_results.append(output)

    plot_data = _plot_study_results(study_results, x=params, y=["Uz", "Ry"])
    print(plot_data)
    print("Parametric study completed.")


def parametric_study_mesh(inputs):
    """Parametric study exploring the effect of material properties on model deflections."""

    study_results = []
    params = [0.5, 1.0, 2.0]  # mesh refinement factor

    for iteration, param_value in enumerate(params):
        print(f"Parametric study iteration {iteration}.")
        inputs_copy = copy.deepcopy(inputs)

        inputs_copy["nele_foil"] = int(
            np.around(inputs_copy["nele_foil"] * param_value)
        )
        # scale the tip loading by the number of airfoil nodes
        new_master = inputs["analysis_file"] + f"_iter{iteration}"
        copy2(inputs["analysis_file"] + ".inp", new_master + ".inp")
        _file_find_replace(
            file=new_master + ".inp",
            find="NLAST,3,1.0",
            replace_with=f"NLAST,3,{1.0/param_value}",
        )
        inputs_copy["analysis_file"] = new_master

        output = main(inputs_copy)
        output["inputs"] = inputs_copy
        study_results.append(output)

    plot_data = _plot_study_results(study_results, x=params, y=["Uz", "Ry"])
    print(plot_data)
    print("Parametric study completed.")


def _plot_study_results(output, x, y):

    y_series = []
    for result in output:
        y_series.append([result[label] for label in y])

    lineObjects = plt.plot(x, y_series)
    plt.xlabel("ply rotation angle (deg)")
    plt.xlabel("displacements (m) or rotation (rad)")
    plt.legend(iter(lineObjects), y)
    plt.show()

    return {
        "x_label": "rotation_angles",
        "x_values": x,
        "y_labels": y,
        "y_values": y_series,
    }


if __name__ == "__main__":
    parametric_study_material(INPUTS[1])
    # parametric_study_mesh(INPUTS[1])
