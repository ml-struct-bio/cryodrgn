import pytest
import os
import argparse
from cryodrgn.source import ImageSource
from cryodrgn.commands import eval_images


@pytest.fixture
def mrcs_data():
    return ImageSource.from_file(f"{pytest.data_dir}/hand.mrcs").images()


@pytest.mark.parametrize(
    "particles, poses, weights, configs",
    [("hand", "hand-rot", "het", "het")],
    indirect=True,
)
def test_invert_contrast(tmpdir, particles, poses, weights, configs):
    args = eval_images.add_args(argparse.ArgumentParser()).parse_args(
        [
            particles.path,
            weights.path,
            "--config",
            configs.path,
            "-o",
            os.path.join(tmpdir, "out_eval_images_losses.pkl"),
            "--out-z",
            os.path.join(tmpdir, "out_eval_images_z.pkl"),
            "--poses",
            poses.path,
            "--log-interval",
            "1",
            "--verbose",
        ]
    )
    eval_images.main(args)
