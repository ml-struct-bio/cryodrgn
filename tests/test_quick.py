import os.path
import argparse
import pytest
from cryodrgn.commands import train_vae


DATA_FOLDER = os.path.join(os.path.dirname(__file__), '..', 'testing', 'data')


@pytest.fixture
def mrcs_file():
    return f'{DATA_FOLDER}/hand.mrcs'


@pytest.fixture
def poses_file():
    return f'{DATA_FOLDER}/hand_rot.pkl'


def test_run(mrcs_file, poses_file):
    args = train_vae.add_args(argparse.ArgumentParser()).parse_args([
        mrcs_file,
        '-o', 'output',
        '--lr', '.0001',
        '--seed', '0',
        '--poses', poses_file,
        '--zdim', '10',
        '--pe-type', 'gaussian'
    ])
    train_vae.main(args)