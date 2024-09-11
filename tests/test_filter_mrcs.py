import pytest
import os
import argparse
import numpy as np
import torch
from cryodrgn.source import ImageSource
from cryodrgn.commands_utils import filter_mrcs
from cryodrgn.utils import save_pkl


@pytest.mark.parametrize("particles", ["toy.mrcs", "hand"], indirect=True)
@pytest.mark.parametrize("ind_size", [0.02, 0.05, 0.25])
@pytest.mark.parametrize("random_seed", [88, 99])
def test_filter_mrcs(tmpdir, particles, ind_size, random_seed):
    mrcs_data = ImageSource.from_file(particles.path).images()

    ind_n = int(mrcs_data.shape[0] * ind_size) + 1
    np.random.seed(random_seed)
    indices = np.random.randint(0, mrcs_data.shape[0], size=ind_n)
    test_lbl = "-".join(
        [particles.label, str(ind_n), str(ind_size)[2:], str(random_seed)]
    )
    ind_fl = os.path.join(tmpdir, f"random-indices_{test_lbl}.pkl")
    save_pkl(indices, ind_fl)

    out_fl = os.path.join(tmpdir, f"projections-filtered_{test_lbl}.mrc")
    args = argparse.ArgumentParser()
    filter_mrcs.add_args(args)
    filter_mrcs.main(args.parse_args([particles.path, "--ind", ind_fl, "-o", out_fl]))

    new_data = ImageSource.from_file(out_fl).images()
    assert torch.allclose(new_data[:], mrcs_data[indices])
