import argparse
import os
import os.path
from cryodrgn.commands import abinit_het, abinit_homo, analyze, backproject_voxel

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")

"""
These tests have unrealistically low parameter values to allow the tests to run on consumer-grade machines quickly,
with or without GPUs.
"""


def test_abinit_het_and_backproject():
    os.makedirs("output", exist_ok=True)

    abinit_args = [
        "--zdim",
        "8",
        "-o",
        "output/abinit_het",
        "--multigpu",  # ok to specify, will be ignored if no GPU(s).
        "--enc-dim",
        "8",
        "--enc-layers",
        "2",
        "--dec-dim",
        "8",
        "--pe-dim",
        "8",
        "--no-window",
    ]

    args = abinit_het.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{DATA_FOLDER}/hand.mrcs",
        ]
        + abinit_args
    )
    abinit_het.main(args)

    # Load a check-pointed model
    abinit_het.main(
        abinit_het.add_args(argparse.ArgumentParser()).parse_args(
            [
                f"{DATA_FOLDER}/hand.mrcs",
                "--load",
                "output/abinit_het/weights.20.pkl",
                "--load-poses",
                "output/abinit_het/pose.20.pkl",
            ]
            + abinit_args
        )
    )

    args = analyze.add_args(argparse.ArgumentParser()).parse_args(
        [
            "output/abinit_het",
            "29",  # Epoch number to analyze - 0-indexed
        ]
    )
    analyze.main(args)

    args = backproject_voxel.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{DATA_FOLDER}/hand.mrcs",
            "--poses",
            "output/abinit_het/pose.pkl",
            "-o",
            "output/abinit_het/backproject/vol.mrc",
        ]
    )
    backproject_voxel.main(args)


def test_abinit_homo_and_backproject():
    os.makedirs("output", exist_ok=True)

    abinit_args = [
        "-o",
        "output/abinit_homo",
        "--dim",
        "16",
        "--pe-dim",
        "8",
    ]

    args = abinit_homo.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{DATA_FOLDER}/hand.mrcs",
        ]
        + abinit_args
    )
    abinit_homo.main(args)

    # Load a check-pointed model
    abinit_homo.main(
        abinit_homo.add_args(argparse.ArgumentParser()).parse_args(
            [
                f"{DATA_FOLDER}/hand.mrcs",
                "--load",
                "output/abinit_homo/weights.20.pkl",
                "--load-poses",
                "output/abinit_homo/pose.20.pkl",
            ]
            + abinit_args
        )
    )

    args = backproject_voxel.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{DATA_FOLDER}/hand.mrcs",
            "--poses",
            "output/abinit_homo/pose.pkl",
            "-o",
            "output/abinit_homo/backproject/vol.mrc",
        ]
    )
    backproject_voxel.main(args)
