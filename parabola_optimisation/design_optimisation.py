"""Design optimisation module making use of OpenMDAO and Scipy libraries. 

Based on the OpenMDAO "Optimization of Paraboloid" example (accessed 09/02/2022)
@ https://openmdao.org/newdocs/versions/latest/basic_user_guide/single_disciplinary_optimization/first_optimization.html """

import openmdao.api as om  # type: ignore
import scipy.optimize.slsqp as scipy_slsqp  # type: ignore

INPUTS = [
    {
        "components": [
            {"name": "Paraboloid", "init": "(inputs,)"},
            {"name": "Constraint", "init": "('g = x + y',)"},
        ],
        "input_data": [
            {
                "component": "Paraboloid",
                "name": "x",
                "value": 50.0,
                "lower": -50,
                "upper": 50,
            },
            {
                "component": "Paraboloid",
                "name": "y",
                "value": 50.0,
                "lower": -50,
                "upper": 50,
            },
        ],
        "connections": [
            {
                "origin": "Paraboloid",
                "name_origin": "x_out",
                "target": "Constraint",
                "name_target": "x",
            },
            {
                "origin": "Paraboloid",
                "name_origin": "y_out",
                "target": "Constraint",
                "name_target": "y",
            },
        ],
        "output_data": [
            {
                "component": "Paraboloid",
                "type": "internal",
                "name": "x_out",
                "value": 0.0,
            },
            {
                "component": "Paraboloid",
                "type": "internal",
                "name": "y_out",
                "value": 0.0,
            },
            {
                "component": "Paraboloid",
                "type": "objective",
                "name": "f_xy",
                "value": 0.0,
            },
            {
                "component": "Constraint",
                "type": "constraint",
                "name": "g",
                "lower": 0.0,
                "upper": 10.0,
            },
        ],
        "optimisation_options": {
            "optimizer": "SLSQP",
            "maxiter": 20,
            "tol": 1e-8,
            "disp": True,
            "debug_print": ["desvars", "ln_cons", "nl_cons", "objs", "totals"],
        },
    }
]


class Paraboloid(om.ExplicitComponent):
    """
    Evaluates the equation f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3.
    """

    def __init__(self, component):
        self.component_inputs = component
        super().__init__()

    def setup(self):

        inputs = self.component_inputs["input_data"]
        outputs = self.component_inputs["output_data"]

        for input in inputs:
            if input["component"] == "Paraboloid":
                self.add_input(input["name"], val=input["value"])

        for output in outputs:
            if output["component"] == "Paraboloid":
                self.add_output(output["name"], val=output["value"])

    def setup_partials(self):
        # Finite difference all partials.
        self.declare_partials("*", "*", method="fd")

    def compute(self, inputs, outputs):
        """
        f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3

        Minimum at: x = 6.6667; y = -7.3333
        """
        x = inputs["x"]
        y = inputs["y"]

        outputs["f_xy"] = (x - 3.0) ** 2 + x * y + (y + 4.0) ** 2 - 3.0
        outputs["x_out"] = x
        outputs["y_out"] = y


component_lookup = {"Paraboloid": Paraboloid, "Constraint": om.ExecComp}


def main(inputs):

    # 1) define the simulation components
    prob = om.Problem()
    for component in inputs["components"]:
        prob.model.add_subsystem(
            component["name"],
            component_lookup[component["name"]](*eval(component["init"])),
        )

    # 2) define the internal component connections
    for connection in inputs["connections"]:
        prob.model.connect(
            connection["origin"] + "." + connection["name_origin"],
            connection["target"] + "." + connection["name_target"],
        )

    # 3) setup the optimisation driver options
    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options["optimizer"] = inputs["optimisation_options"]["optimizer"]
    prob.driver.options["maxiter"] = inputs["optimisation_options"]["maxiter"]
    prob.driver.options["tol"] = inputs["optimisation_options"]["tol"]
    prob.driver.options["disp"] = inputs["optimisation_options"]["disp"]
    prob.driver.options["debug_print"] = inputs["optimisation_options"]["debug_print"]
    # Ask OpenMDAO to finite-difference across the model to compute the gradients for the optimizer
    prob.model.approx_totals()  # this forces FD gradients

    # 4) add design variables
    for var in inputs["input_data"]:
        prob.model.add_design_var(
            f"{var['component']}.{var['name']}", lower=var["lower"], upper=var["upper"]
        )

    # 5) add an objective and constraints
    for var in inputs["output_data"]:
        if var["type"] == "objective":
            prob.model.add_objective(f"{var['component']}.{var['name']}")
        elif var["type"] == "constraint":
            prob.model.add_constraint(
                f"{var['component']}.{var['name']}",
                lower=var["lower"],
                upper=var["upper"],
            )

    # 6) execute the optimisation
    prob.setup()
    om.n2(prob, outfile="n2.html")  # visualise the n2 diagram
    print("=== optimisation started ===")
    prob.run_driver()
    print("=== optimisation stopped ===")

    # print objective and constraint values
    for var in inputs["output_data"]:
        if var["type"] == "objective":
            name = f"{var['component']}.{var['name']}"
            print("Objective " + name + " value: ", prob.get_val(name))
        if var["type"] == "constraint":
            name = f"{var['component']}.{var['name']}"
            print("Constraint " + name + " value: ", prob.get_val(name))

    # print variable values
    for var in inputs["input_data"]:
        name = f"{var['component']}.{var['name']}"
        print("Variable " + name + " value: ", prob.get_val(name))


if __name__ == "__main__":
    main(INPUTS[0])