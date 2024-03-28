import os.path
from cryodrgn.utils import run_command

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


def test_fidelity_small():
    zvals_fl = os.path.join(DATA_FOLDER, "zvals_het-2_1k.pkl")
    out, err = run_command(
        f"cryodrgn graph_traversal {zvals_fl} --anchors 50 100 "
        f"--max-neighbors=30 --avg-neighbors=10"
    )
    assert err == ""

    assert "\n".join(out.split("\n")[1:6] + out.split("\n")[7:]) == (
        "Working on images 0-1000\n"
        "Working on images 1000-1319\n"
        "Max dist between neighbors: 0.1022  (to enforce average of 10.0 neighbors)\n"
        "Total path distance 50->100: 0.8441\n"
        " Euclidean distance 50->100: 0.7868\n"
        "ind	z1	z2	dist\n"
        "50	-1.974	-0.253	0.0\n"
        "22	-1.908	-0.18	0.099\n"
        "1310	-1.838	-0.149	0.076\n"
        "1047	-1.745	-0.124	0.096\n"
        "1093	-1.673	-0.145	0.075\n"
        "857	-1.59	-0.119	0.087\n"
        "159	-1.488	-0.114	0.102\n"
        "1215	-1.44	-0.074	0.063\n"
        "157	-1.378	-0.068	0.062\n"
        "1006	-1.325	-0.083	0.055\n"
        "672	-1.265	-0.089	0.06\n"
        "100	-1.216	-0.04	0.069\n\n"
    )


def test_fidelity_medium():
    zvals_fl = os.path.join(DATA_FOLDER, "zvals_het-2_1k.pkl")
    out, err = run_command(
        f"cryodrgn graph_traversal {zvals_fl} --anchors 50 100 "
        f"--max-neighbors=50 --avg-neighbors=20"
    )
    assert err == ""

    assert "\n".join(out.split("\n")[1:6] + out.split("\n")[7:]) == (
        "Working on images 0-1000\n"
        "Working on images 1000-1319\n"
        "Max dist between neighbors: 0.149  (to enforce average of 20.0 neighbors)\n"
        "Total path distance 50->100: 0.808\n"
        " Euclidean distance 50->100: 0.7868\n"
        "ind	z1	z2	dist\n"
        "50	-1.974	-0.253	0.0\n"
        "22	-1.908	-0.18	0.099\n"
        "1310	-1.838	-0.149	0.076\n"
        "1047	-1.745	-0.124	0.096\n"
        "1069	-1.638	-0.093	0.111\n"
        "688	-1.492	-0.084	0.147\n"
        "913	-1.35	-0.061	0.144\n"
        "100	-1.216	-0.04	0.135\n\n"
    )
