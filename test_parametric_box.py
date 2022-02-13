""" Test functions for parametric_box module."""

import os
import sys
import pytest
import numpy as np
from pathlib import Path

from parametric_box import _get_from_dat, _get_average_rotation

_TEST_DATA_FOLDER = Path("test_data")


@pytest.mark.parametrize(
    "cases",
    [
        {
            "output_file": "ccx_static_tip_shear_Rx.dat",
            "mesh_file": "all.msh",
            "Rx": 0.01,
        },
        {
            "output_file": "ccx_static_tip_shear_Ry.dat",
            "mesh_file": "all.msh",
            "Ry": 0.01,
        },
        {
            "output_file": "ccx_static_tip_shear_Rz.dat",
            "mesh_file": "all.msh",
            "Rz": 0.01,
        },
    ],
)
def test_get_average_rotation(cases):

    # read file and recover node displacements
    output_file = str(_TEST_DATA_FOLDER / "rigid_body_rotation" / cases["output_file"])
    # all_disp : [nodeid, vx, vy, vz]
    all_disps = _get_from_dat(output_file, data="displacements")

    # calculate the average rotations
    mesh_file = str(_TEST_DATA_FOLDER / "rigid_body_rotation" / cases["mesh_file"])
    rotations_mean = _get_average_rotation(all_disps=all_disps, mesh_file=mesh_file)

    if "Rx" in cases:
        assert np.isclose(cases["Rx"], rotations_mean[0], rtol=1e-3)
    if "Ry" in cases:
        assert np.isclose(cases["Ry"], rotations_mean[1], rtol=1e-3)
    if "Rz" in cases:
        assert np.isclose(cases["Rz"], rotations_mean[2], rtol=1e-3)
