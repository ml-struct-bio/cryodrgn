import numpy as np
from numpy.testing import assert_array_almost_equal

from cryodrgn import utils


def test_convert_from_relion_scipy():
    x = np.array([[-147.9485, 22.58999, 162.539]])
    r1 = utils.R_from_relion_scipy(x)[0]
    r2 = utils.R_from_relion(*x[0])
    assert_array_almost_equal(r1, r2)


def test_convert_from_relion():
    x = np.array([[-147.9485, 22.58999, 162.539]])
    y = np.array(
        [
            [
                [0.90571941, 0.21306972, 0.36643367],
                [-0.27142096, 0.95553406, 0.11526193],
                [-0.32558103, -0.20385275, 0.92327734],
            ]
        ]
    )
    r1 = utils.R_from_relion_scipy(x)
    assert_array_almost_equal(r1, y)


def test_convert_to_relion():
    x = np.array([[-147.9485, 22.58999, 162.539]])
    r1 = utils.R_from_relion_scipy(x)
    euler = utils.R_to_relion_scipy(r1)
    assert_array_almost_equal(x, euler)
