import pytest
import os.path
from cryodrgn.utils import run_command


def test_fidelity_small():
    zvals_fl = os.path.join(pytest.DATADIR, "zvals_het-2_1k.pkl")
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
        "50(a)	-1.974	-0.253	0.0\n"
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
        "100(a)	-1.216	-0.04	0.069\n\n"
    )


def test_fidelity_medium():
    zvals_fl = os.path.join(pytest.DATADIR, "zvals_het-2_1k.pkl")
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
        "50(a)	-1.974	-0.253	0.0\n"
        "22	-1.908	-0.18	0.099\n"
        "1310	-1.838	-0.149	0.076\n"
        "1047	-1.745	-0.124	0.096\n"
        "1069	-1.638	-0.093	0.111\n"
        "688	-1.492	-0.084	0.147\n"
        "913	-1.35	-0.061	0.144\n"
        "100(a)	-1.216	-0.04	0.135\n\n"
    )


def test_fidelity_large():
    zvals_fl = os.path.join(pytest.DATADIR, "zvals_het-8_4k.pkl")
    out, err = run_command(
        f"cryodrgn graph_traversal {zvals_fl} --anchors 2020 3030 4040 "
        f"--max-neighbors=10 --avg-neighbors=5"
    )
    assert err == ""

    assert "\n".join(out.split("\n")[1:11] + out.split("\n")[12:]) == (
        "Working on images 0-1000\n"
        "Working on images 1000-2000\n"
        "Working on images 2000-3000\n"
        "Working on images 3000-4000\n"
        "Working on images 4000-4549\n"
        "Max dist between neighbors: 0.7912  (to enforce average of 5.0 neighbors)\n"
        "Total path distance 2020->3030: 2.409\n"
        " Euclidean distance 2020->3030: 1.528\n"
        "Total path distance 3030->4040: 5.141\n"
        " Euclidean distance 3030->4040: 2.689\n"
        "ind	z1	z2	z3	z4	z5	z6	z7	z8	dist\n"
        "2020(a)	-2.132	-0.063	1.11	-0.422	0.539	-0.169	0.187	0.2	0.0\n"
        "4117	-2.108	-0.115	0.939	-0.465	0.325	-0.154	0.233	0.362	0.33\n"
        "1402	-2.112	-0.019	0.621	-0.615	0.335	-0.056	0.212	0.313	0.381\n"
        "3480	-2.397	-0.137	0.157	-0.524	0.249	0.077	0.217	0.397	0.592\n"
        "1106	-2.582	-0.053	0.102	-0.193	0.13	-0.047	0.112	0.787	0.589\n"
        "3030(a)	-2.616	-0.06	0.14	-0.064	0.302	-0.009	-0.263	1.064	0.517\n"
        "3824	-2.964	0.082	0.12	0.006	0.236	-0.04	0.058	0.897	0.532\n"
        "3792	-2.583	0.124	0.149	0.28	0.277	-0.021	-0.056	0.297	0.773\n"
        "707	-2.553	0.09	-0.191	0.354	0.394	0.041	-0.06	0.151	0.402\n"
        "41	-2.401	0.274	-0.211	0.807	0.234	0.043	0.219	0.074	0.609\n"
        "3210	-2.249	0.17	-0.605	0.652	-0.1	0.091	0.0	0.073	0.612\n"
        "2257	-2.236	-0.028	-0.886	0.829	-0.059	0.198	0.012	0.218	0.43\n"
        "392	-2.491	0.09	-1.328	0.768	0.05	-0.031	-0.019	0.383	0.609\n"
        "301	-2.669	0.139	-1.748	0.938	-0.215	-0.042	-0.187	0.528	0.599\n"
        "4040(a)	-2.531	0.15	-2.177	0.769	-0.414	-0.053	-0.064	0.315	0.575\n\n"
    )
