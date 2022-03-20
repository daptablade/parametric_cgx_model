""" Static strip theory aerodynamics for symmetric aerofoils implemented in python."""

from re import M
import numpy as np
import matplotlib.pyplot as plt

from parametric_box import INPUTS, PLOT_FLAG

AERO_INPUTS = {
    "planform": None,  # from box model - or {"A":[0,0,0], "B":[0,2.0,0], "CA":0.2, "CB":0.2}
    "strips": 10,  # number of aero strips along the span >=1
    "root_alpha": 10,  # AoA at the wing root in degrees
    "rho": 1.225,  # air density in Pa
    "V": 10,  # air velocity in m/s
    "CL_alpha": 2 * np.pi,  # lift curve slope
}


def plot_forces(nodes, forces, le_nodes, te_nodes):
    """Plot the force distribution in 3D over the aero planform."""

    max_force = np.max(np.linalg.norm(forces, axis=1))
    nforces = forces / max_force * 0.1
    ax = plt.figure().add_subplot(projection="3d")
    for (le, te) in zip(le_nodes, te_nodes):
        ax.plot([le[0], te[0]], [le[1], te[1]], [le[2], te[2]], "o-b")
    ax.quiver3D(
        nodes[:, 0],
        nodes[:, 1],
        nodes[:, 2],
        nforces[:, 0],
        nforces[:, 1],
        nforces[:, 2],
        normalize=False,
    )
    plt.show()


class AeroModel(object):
    """Aerodynamic model class."""

    def __init__(self, aero_inputs, box_inputs):
        self.aero_inputs = aero_inputs
        self.box_inputs = box_inputs
        self._point_a = None
        self._point_b = None
        self._chord_a = None
        self._chord_b = None
        self.le_nodes = None
        self.te_nodes = None  # same length as le_nodes
        self.alpha = None  # same length as le_nodes
        self.lift = None  # length of le_nodes -1
        self.le_aeroforces = None  # same length as le_nodes
        self.te_aeroforces = None  # same length as le_nodes

    def get_planform(self):
        """The planform is defined by points at the leading edge and root (A),
        and leading edge and tip (B) of the wing and the aero chords at the root (CA),
        and a the tip (CB).
        """

        if not self.aero_inputs["planform"]:
            span = self.box_inputs["span"]
            chord = self.box_inputs["chord"]
            self._point_a = np.array([0.0, 0.0, 0.0], dtype=float)
            self._point_b = np.array([0.0, span, 0.0], dtype=float)
            self._chord_a = chord
            self._chord_b = chord
        else:
            self._point_a = np.array(self.aero_inputs["planform"]["A"], dtype=float)
            self._point_b = np.array(self.aero_inputs["planform"]["B"], dtype=float)
            self._chord_a = self.aero_inputs["planform"]["CA"]
            self._chord_b = self.aero_inputs["planform"]["CB"]

    def get_all_nodes(self):
        """calculate the internal node positions at the leading and trailing edges."""

        linear_int = lambda a, b, n, ii: (b - a) / n * ii + a
        strips = self.aero_inputs["strips"]
        self.le_nodes = np.zeros([strips + 1, 3])
        self.te_nodes = np.zeros([strips + 1, 3])
        for inc in range(strips + 1):
            for coord in range(3):
                self.le_nodes[inc, coord] = linear_int(
                    self._point_a[coord], self._point_b[coord], strips, inc
                )
            chord = linear_int(self._chord_a, self._chord_b, strips, inc)
            self.te_nodes[inc, :] = self.le_nodes[inc, :] + np.array([chord, 0.0, 0.0])

    def get_alphas(self, node_displacements=None):
        """Calculate the effective angle of incidence at each strip edge."""

        alpha_flex = np.zeros(
            [len(self.le_nodes), 1]
        )  # TODO calculate from node_displacements

        self.alpha = alpha_flex + self.aero_inputs["root_alpha"]

    def compute_lift(self):
        """Compute the lift force at each aero strip."""

        lift = []
        strips = self.aero_inputs["strips"]
        rho = self.aero_inputs["rho"]
        v = self.aero_inputs["V"]
        cl_a = self.aero_inputs["CL_alpha"]

        for strip in range(strips):
            a_root = self.alpha[strip]
            a_tip = self.alpha[strip + 1]
            alpha_mean = (a_root + a_tip) / 2

            c_root = self.te_nodes[strip, 0] - self.le_nodes[strip, 0]
            c_tip = self.te_nodes[strip + 1, 0] - self.le_nodes[strip + 1, 0]
            delta_span = self.le_nodes[strip + 1, 1] - self.le_nodes[strip, 1]
            ds = (c_root + c_tip) / 2 * delta_span
            strip_lift = 0.5 * rho * v**2 * cl_a * alpha_mean * ds
            lift.append(strip_lift)

        self.lift = np.array(lift, dtype=float)

    def compute_aeroforces_at_nodes(self):
        """Compute the aero forces at each mesh node."""

        strips = self.aero_inputs["strips"]
        lift = np.zeros([strips + 1, 3])
        for strip in range(strips):
            lift[strip, 2] += self.lift[strip] / 2
            lift[strip + 1, 2] += self.lift[strip] / 2

        # aerodynamic center assumed to be at 1/4 chord
        self.le_aeroforces = lift * 3 / 4
        self.te_aeroforces = lift * 1 / 4


def main(box_inputs):
    """Create the aero model and calculate the nodal forces as a function of
    the flight conditions.

    rho = air density
    V = air velocity
    alpha = effective angle of attack
    CL_alpha = lift curve slope
    """

    # instantiate the aeromodel class
    aeromodel = AeroModel(AERO_INPUTS, box_inputs)

    # define the corner nodes of the aerodynamic planform
    aeromodel.get_planform()

    # calculate the internal node positions at the leading and trailing edges
    aeromodel.get_all_nodes()

    # recover the local average angle of attack for each strip
    aeromodel.get_alphas()

    # calculate the Lift at the aerodynamic center of each strip
    aeromodel.compute_lift()

    # calculate the equivalent forces at the leading and trailing edges of the strip
    aeromodel.compute_aeroforces_at_nodes()

    aeronodes = np.vstack([aeromodel.le_nodes, aeromodel.te_nodes])
    aeroforces = np.vstack([aeromodel.le_aeroforces, aeromodel.te_aeroforces])

    # plot the force distribution at the aero nodes
    if PLOT_FLAG:
        plot_forces(
            nodes=aeronodes,
            forces=aeroforces,
            le_nodes=aeromodel.le_nodes,
            te_nodes=aeromodel.te_nodes,
        )

    return aeronodes, aeroforces


if __name__ == "__main__":
    _ = main(INPUTS[1])
