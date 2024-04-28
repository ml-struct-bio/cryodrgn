import os
import numpy as np
from cryodrgn.utils import run_command

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


def test_fidelity_small():
    out_file = os.path.join("output", "path.txt")
    out, err = run_command(
        f"cryodrgn direct_traversal {os.path.join(DATA_FOLDER, 'zvals_het-2_1k.pkl')} "
        f"--anchors 77 333 -o {out_file}"
    )
    assert err == ""

    path = np.loadtxt(out_file)
    assert np.allclose(
        path,
        np.array(
            [
                [-2.08011341, -0.1485731],
                [-1.96137881, -0.14291579],
                [-1.84264421, -0.13725847],
                [-1.72390962, -0.13160115],
                [-1.60517502, -0.12594384],
                [-1.48644042, -0.12028652],
            ]
        ),
    )


def test_fidelity_big():
    out_file = os.path.join("output", "path.txt")
    out, err = run_command(
        f"cryodrgn direct_traversal {os.path.join(DATA_FOLDER, 'zvals_het-8_4k.pkl')} "
        f"--anchors 1099 4001 -o {out_file}"
    )
    assert err == ""

    path = np.loadtxt(out_file)
    assert round((path @ path.T).sum(), 5) == 185.24840
