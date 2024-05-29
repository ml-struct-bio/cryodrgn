import pytest
import os.path
from cryodrgn.utils import run_command


def test_fidelity_small():
    out, err = run_command(
        f"cryodrgn pc_traversal {os.path.join(pytest.DATADIR, 'zvals_het-2_1k.pkl')} "
        "--pc 0 --lim 0.10 0.85 -n 5"
    )
    assert err == ""

    outs = out.split("\n")
    limit_txts = outs[3].split(": ")[-1].split(", ")
    assert round(float(limit_txts[0]), 5) == -0.83899
    assert round(float(limit_txts[1]), 5) == -0.64078
    assert outs[5][1:-1].split() == ["1183", "1186", "1189", "1191", "1198"]


def test_fidelity_big():
    out, err = run_command(
        f"cryodrgn pc_traversal "
        f"{os.path.join(pytest.DATADIR, 'zvals_het-8_4k.pkl')} --pc 3"
    )
    assert err == ""

    outs = out.split("\n")
    limit_txts = outs[4].split(": ")[-1].split(", ")
    assert round(float(limit_txts[0]), 5) == -2.17560
    assert round(float(limit_txts[1]), 5) == 2.29041

    assert outs[6][1:-1].split() == [
        "1016",
        "1324",
        "1611",
        "1823",
        "1951",
        "1905",
        "1686",
        "1407",
        "1043",
        "673",
    ]
