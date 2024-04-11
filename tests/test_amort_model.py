"""Fidelity tests of the hierarchical pose search reconstruction algorithm.

"""
import pytest
import os
from cryodrgn.utils import run_command

REFVOL_DIR = os.path.join(pytest.data_dir, "reference-volumes")


def same_volumes(vol1, vol2):
    fsc_out, err = run_command(f"cryodrgn_utils fsc {vol1} {vol2}")

    fsc_val = 0.0
    for fsc_str in fsc_out.split("\n"):
        if fsc_str[:2] == "0." and len(fsc_str.split()) > 1:
            fsc_val = float(fsc_str.split()[1])

    return fsc_val > 0.999


@pytest.mark.parametrize(
    "reconstruct_dir",
    [
        {
            "dataset": "hand_16",
            "model": "amort",
            "seed": 8543,
            "hidden_dim": 24,
            "pe_dim": 8,
            "epochs_sgd": 2,
            "n_imgs_pretrain": 100,
            "n_imgs_pose_search": 200,
            "t_extent": 4.0,
            "t_ngrid": 2,
            "num_workers": 1,
            "amp": False,
        }
    ],
    indirect=True,
)
def test_homo_abinit(reconstruct_dir):
    reconstruct_dir.train()
    out_mrc = os.path.join(reconstruct_dir.outdir, "out", "reconstruct.5.mrc")
    assert os.path.exists(out_mrc)

    comp_mrc = os.path.join(REFVOL_DIR, "abinit-homo_simple.mrc")
    assert same_volumes(comp_mrc, out_mrc)
