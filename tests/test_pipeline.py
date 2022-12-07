import argparse
import os
import os.path
import pytest
from cryodrgn.commands import downsample, parse_ctf_star, parse_pose_star
from cryodrgn.commands_utils import write_star
from cryodrgn.utils import assert_pkl_close

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def toy_projections_star():
    return f"{DATA_FOLDER}/toy_projections.star"


@pytest.fixture
def toy_projections_txt():
    return f"{DATA_FOLDER}/toy_projections.txt"


@pytest.fixture
def toy_particles_mrcs():
    return f"{DATA_FOLDER}/toy_projections.mrcs"


@pytest.fixture
def relion_mrcs():
    return f"{DATA_FOLDER}/relion31.mrcs"


@pytest.fixture
def relion_starfile():
    return f"{DATA_FOLDER}/FinalRefinement-OriginalParticles-PfCRT.star"


def test_pipeline(
    toy_projections_star,
    toy_projections_txt,
    toy_particles_mrcs,
    relion_mrcs,
    relion_starfile,
):
    os.makedirs("output", exist_ok=True)

    # Downsample all particles from an .mrcs file to a pkl file
    args = parse_ctf_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            toy_projections_star,
            "-D",
            "30",
            "--Apix",
            "1",
            "-o",
            "output/pipeline_ctf.pkl",
        ]
    )
    parse_ctf_star.main(args)

    # Write starfile from an input .mrcs, with ALL particles selected
    args = write_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            toy_particles_mrcs,
            "--ctf",
            "output/pipeline_ctf.pkl",
            "-o",
            "output/pipeline_test.star",
        ]
    )
    write_star.main(args)

    # Write starfile from an input .mrcs, with 100 particles selected
    args = write_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            toy_particles_mrcs,
            "--ctf",
            "output/pipeline_ctf.pkl",
            "--ind",
            f"{DATA_FOLDER}/ind100.pkl",
            "-o",
            "output/pipeline_test100.star",
        ]
    )
    write_star.main(args)

    # Write starfile from an input .mrcs, with poses, with ALL particles selected
    args = write_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            toy_particles_mrcs,
            "--ctf",
            "output/pipeline_ctf.pkl",
            "--poses",
            f"{DATA_FOLDER}/toy_rot_trans.pkl",
            "-o",
            "output/pipeline_test.star",
        ]
    )
    write_star.main(args)

    # Write poses pkl from starfile with ALL particles selected
    args = parse_pose_star.add_args(argparse.ArgumentParser()).parse_args(
        ["output/pipeline_test.star", "-D", "30", "-o", "output/test_pose.pkl"]
    )
    parse_pose_star.main(args)
    assert_pkl_close("output/test_pose.pkl", f"{DATA_FOLDER}/toy_rot_trans.pkl")

    # Downsample all particles from an starfile to an .mrcs file
    args = downsample.add_args(argparse.ArgumentParser()).parse_args(
        [
            toy_projections_star,
            "-D",
            "28",
            "--chunk",
            "80",
            "-o",
            "output/toy_projections.mrcs",
        ]
    )
    downsample.main(args)

    # Write starfile from an input .txt file with 100 particles selected
    args = write_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            toy_projections_txt,
            "--ctf",
            "output/pipeline_ctf.pkl",
            "--ind",
            f"{DATA_FOLDER}/ind100.pkl",
            "-o",
            "output/pipeline_test2.star",
        ]
    )
    write_star.main(args)

    # Test copying micrograph coordinates
    args = write_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            relion_mrcs,
            "--ctf",
            f"{DATA_FOLDER}/ctf1.pkl",
            "-o",
            "output/pipeline_test3.star",
        ]
    )
    write_star.main(args)

    # Test copying micrograph coordinates, with filtering
    args = write_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            relion_mrcs,
            "--ctf",
            f"{DATA_FOLDER}/ctf1.pkl",
            "--ind",
            f"{DATA_FOLDER}/ind3.pkl",
            "-o",
            "output/pipeline_test3_filtered.star",
        ]
    )
    write_star.main(args)
