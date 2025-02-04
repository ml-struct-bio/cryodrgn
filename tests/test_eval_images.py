import pytest
import os
import argparse
from cryodrgn.commands import eval_images


@pytest.mark.parametrize(
    "particles, poses, weights, configs",
    [("hand", "hand-rot", "het", "het"), ("hand", "hand-poses", "het", "het")],
    indirect=True,
)
def test_invert_contrast(tmpdir, particles, poses, weights, configs):
    parser = argparse.ArgumentParser()
    eval_images.add_args(parser)
    args = parser.parse_args(
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
