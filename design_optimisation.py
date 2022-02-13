"""Design optimisation module making use of OpenMDAO and Scipy libraries. """

import numpy as np
from matplotlib import pyplot as plt  # type: ignore

import openmdao.api as om  # type: ignore
import scipy.optimize.slsqp as scipy_slsqp  # type: ignore

import parametric_box as parabox

INPUTS = [
    {
        "components": [
            {"name": "ParaBox", "init": "inputs, parabox.INPUTS[1]"},
        ],
        "input_data": [
            {
                "component": "ParaBox",
                "name": "angle",
                "value": -20.0,
                "lower": -50,
                "upper": 50,
            }
        ],
        "connections": [],
        "output_data": [
            {
                "component": "ParaBox",
                "type": "objective",
                "name": "Ry",
                "value": 0.0,
            },
            {
                "component": "ParaBox",
                "type": "constraint",
                "name": "Uz",
                "value": 0.0,
                "lower": -0.1,
                "upper": 0.1,
            },
        ],
        "optimisation_options": {
            "optimizer": "SLSQP",
            "maxiter": 20,
            "tol": 1e-9,
            "disp": True,
            "debug_print": ["desvars", "ln_cons", "nl_cons", "objs", "totals"],
        },
        "options": {"show_n2_diagram": False},
    }
]


class ParaBox(om.ExplicitComponent):
    """
    Evaluates the static deflections of the box model.
    """

    def __init__(self, component, parameters):
        self.component_inputs = component
        self.component_parameters = parameters
        super().__init__()

    def setup(self):

        inputs = self.component_inputs["input_data"]
        outputs = self.component_inputs["output_data"]

        for input in inputs:
            if input["component"] == "ParaBox":
                self.add_input(input["name"], val=input["value"])

        for output in outputs:
            if output["component"] == "ParaBox":
                self.add_output(output["name"], val=output["value"])

    def setup_partials(self):
        # Finite difference all partials.
        self.declare_partials("*", "*", method="fd")

    def compute(self, inputs, outputs):

        # define optimisation inputs
        angle = inputs["angle"]

        # perform computations of output functions
        def compute_outputs(params, angle):
            params["orientations"][0]["1"] = parabox._rotate_vector(
                angle=angle,
                starting=[0.0, 1.0, 0.0],
                axis=[0.0, 0.0, 1.0],
            )
            return parabox.main(params)

        computed_data = compute_outputs(self.component_parameters, angle)

        # recover optimisation outputs
        outputs["Uz"] = computed_data["Uz"]
        outputs["Ry"] = computed_data["Ry"]


component_lookup = {"ParaBox": ParaBox}


def main(inputs):

    # 1) define the simulation components
    prob = om.Problem()
    for component in inputs["components"]:
        prob.model.add_subsystem(
            component["name"],
            component_lookup[component["name"]](*eval(component["init"])),
        )

    # 2) define the internal component connections => none here
    for connection in inputs["connections"]:
        prob.model.connect(
            connection["origin"] + "." + connection["name_origin"],
            connection["target"] + "." + connection["name_target"],
        )

    # 3) setup the optimisation driver options
    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options["optimizer"] = inputs["optimisation_options"]["optimizer"]
    prob.driver.opt_settings["maxiter"] = inputs["optimisation_options"]["maxiter"]
    prob.driver.opt_settings["ftol"] = inputs["optimisation_options"]["tol"]
    prob.driver.opt_settings["disp"] = inputs["optimisation_options"]["disp"]
    prob.driver.options["debug_print"] = inputs["optimisation_options"]["debug_print"]

    # add design variables
    for var in inputs["input_data"]:
        prob.model.add_design_var(
            f"{var['component']}.{var['name']}", lower=var["lower"], upper=var["upper"]
        )

    # add an objective and constraints
    for var in inputs["output_data"]:
        if var["type"] == "objective":
            prob.model.add_objective(f"{var['component']}.{var['name']}")
        elif var["type"] == "constraint":
            prob.model.add_constraint(
                f"{var['component']}.{var['name']}",
                lower=var["lower"],
                upper=var["upper"],
            )

    # Ask OpenMDAO to finite-difference across the model to compute the gradients for the optimizer
    prob.model.approx_totals(
        method="fd", step=0.1, form="forward", step_calc="abs"
    )  # this forces FD gradients

    # add a data recorder to the optimisation problem
    r = om.SqliteRecorder("problem_recorder.sqlite")
    prob.driver.add_recorder(r)
    prob.driver.recording_options["record_derivatives"] = True

    # execute the optimisation
    prob.setup()
    if inputs["options"]["show_n2_diagram"]:
        om.n2(prob, outfile="n2.html")  # visualise the n2 diagram
    print("=== optimisation started ===")
    prob.run_driver()
    print("=== optimisation stopped ===")

    # print final objective and constraint values
    for var in inputs["output_data"]:
        if var["type"] == "objective":
            name = f"{var['component']}.{var['name']}"
            print("Objective " + name + " value: ", prob.get_val(name))
        if var["type"] == "constraint":
            name = f"{var['component']}.{var['name']}"
            print("Constraint " + name + " value: ", prob.get_val(name))

    # print final variable values
    for var in inputs["input_data"]:
        name = f"{var['component']}.{var['name']}"
        print("Variable " + name + " value: ", prob.get_val(name))


def post_process(inputs, only_plot_major_iter=True):
    # read database
    # Instantiate your CaseReader
    cr = om.CaseReader("problem_recorder.sqlite")
    # Isolate "problem" as your source
    driver_cases = cr.list_cases("driver", out_stream=None)

    # plot the iteration history from the recorder data
    inputs_history = {
        key: []
        for key in [f"{var['component']}.{var['name']}" for var in inputs["input_data"]]
    }
    outputs_history = {
        key: []
        for key in [
            f"{var['component']}.{var['name']}" for var in inputs["output_data"]
        ]
    }
    for key in driver_cases:
        case = cr.get_case(key)
        if (only_plot_major_iter and case.derivatives) or not only_plot_major_iter:
            # get history of inputs
            for key in inputs_history:
                inputs_history[key].append(case.outputs[key])
            # get history of outputs
            for key in outputs_history:
                outputs_history[key].append(case.outputs[key])

    # plot output in userfriendly fashion
    _plot_iteration_histories(
        inputs_history=inputs_history, outputs_history=outputs_history
    )


def _plot_iteration_histories(inputs_history=None, outputs_history=None):
    # plot input histories
    for key in inputs_history:
        input_data = inputs_history[key]
        input_data = np.array(input_data)
        iterations = range(input_data.shape[0])

        plt.figure()
        for data_series in input_data.T:
            plt.plot(iterations, data_series, "-o")
            plt.grid(True)
            plt.title(key)

    # plot output histories
    for key in outputs_history:
        output_data = outputs_history[key]
        output_data = np.array(output_data)
        iterations = range(output_data.shape[0])

        plt.figure()
        for data_series in output_data.T:
            plt.plot(iterations, data_series, "-o")
            plt.grid(True)
            plt.title(key)

    plt.show()


if __name__ == "__main__":
    main(INPUTS[0])
    post_process(INPUTS[0])