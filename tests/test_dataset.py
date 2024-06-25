import pytest
import os
import numpy as np
from typing import Sequence
from torch.utils.data.sampler import BatchSampler, RandomSampler
from torch.utils.data import DataLoader
from cryodrgn.dataset import DataShuffler, ImageDataset, TiltSeriesData, make_dataloader
from cryodrgn.utils import load_pkl


@pytest.mark.parametrize("particles", ["hand"], indirect=True)
def test_particles(particles):
    # compare data for some random indices against data generated from known
    # correct values after normalization/fft/symmetrization etc.
    dataset = ImageDataset(mrcfile=particles.path, invert_data=True)
    indices = np.array([43, 12, 53, 64, 31, 56, 75, 63, 27, 62, 96])
    particles_arr, _, _ = dataset[indices]

    # the tolerance (atol) of 1e-5 accounts for minor fft differences between
    # numpy fft (legacy) and torch fft
    assert np.allclose(
        np.load(os.path.join(pytest.DATADIR, "hand_11_particles.npy")),
        particles_arr.cpu().numpy(),
        atol=1e-5,
    )


@pytest.mark.parametrize(
    "particles, indices",
    [
        ("hand", None),
        ("hand-tilt", None),
        ("toy.mrcs", None),
        ("toy.mrcs", "first-100"),
        ("toy.mrcs", "random-100"),
        ("toy.txt", None),
        ("toy.txt", "random-100"),
        ("toy.star", None),
        ("toy.star", "random-100"),
        ("tilts.star", None),
    ],
    indirect=True,
)
class TestImageDatasetLoading:

    # 11 is a useful number to test as it is 1 mod 100 and 1 mod 1000
    # lots of edge cases with numpy and torch arrays when a dimension is of length one!
    @pytest.mark.parametrize("batch_size", [7, 11, 20, 43])
    def test_loading_slow(self, particles, indices, batch_size):
        ind = None if indices.path is None else load_pkl(indices.path)
        dataset = ImageDataset(mrcfile=particles.path, ind=ind)
        data_loader = make_dataloader(dataset, batch_size=batch_size, shuffle=True)

        # minibatch is a list of (particles, tilt, indices)
        for i, minibatch in enumerate(data_loader):
            assert isinstance(minibatch, Sequence)
            assert len(minibatch) == 3

            # We have 100 particles. For all but the last iteration *
            # for all but the last iteration (100//7 = 14), we'll have 7 images each
            if i < (dataset.N // batch_size):
                assert minibatch[0].shape == (batch_size, dataset.D, dataset.D)
                assert minibatch[1] is None
                assert minibatch[2].shape == (batch_size,)

            # and 100 % 7 = 2 for the very last one
            else:
                assert minibatch[0].shape == (
                    dataset.N % batch_size,
                    dataset.D,
                    dataset.D,
                )
                assert minibatch[1] is None
                assert minibatch[2].shape == (dataset.N % batch_size,)

    @pytest.mark.parametrize("batch_size", [11, 25, 61])
    def test_loading_fast(self, particles, indices, batch_size):
        ind = None if indices.path is None else load_pkl(indices.path)
        dataset = ImageDataset(mrcfile=particles.path, ind=ind)

        # A faster way to load is to use BatchSampler with RandomSampler
        # see https://stackoverflow.com/questions/61458305
        data_loader = DataLoader(
            dataset,
            sampler=BatchSampler(
                RandomSampler(dataset),
                batch_size=batch_size,
                drop_last=False,  # Or SequentialSampler
            ),
            batch_size=None,  # Do not let torch add a new batch dimension at the front
        )

        # minibatch is a list of (particles, tilt, indices)
        for i, minibatch in enumerate(data_loader):
            assert isinstance(minibatch, Sequence)
            assert len(minibatch) == 3

            # We have 100 particles. For all but the last iteration *
            # for all but the last iteration (100//7 = 14), we'll have 7 images each
            if i < (dataset.N // batch_size):
                assert minibatch[0].shape == (batch_size, dataset.D, dataset.D)
                assert minibatch[1] is None
                assert minibatch[2].shape == (batch_size,)

            # and 100 % 7 = 2 for the very last one
            else:
                assert minibatch[0].shape == (
                    dataset.N % batch_size,
                    dataset.D,
                    dataset.D,
                )

                assert minibatch[1] is None
                assert minibatch[2].shape == (dataset.N % batch_size,)


@pytest.mark.parametrize(
    "particles, indices",
    [
        ("tilts.star", None),
        ("tilts.star", "just-5"),
    ],
    indirect=True,
)
@pytest.mark.parametrize("ntilts", [2, 5, 10])
class TestTiltSeriesLoading:
    @pytest.mark.parametrize("batch_size", [3, 5, 8])
    def test_loading_slow(self, particles, indices, ntilts, batch_size):
        pt_ind = None if indices.path is None else load_pkl(indices.path)

        if pt_ind is None:
            ind = None
        else:
            pt, tp = TiltSeriesData.parse_particle_tilt(particles.path)
            ind = TiltSeriesData.particles_to_tilts(pt, pt_ind)

        dataset = TiltSeriesData(tiltstar=particles.path, ntilts=ntilts, ind=ind)
        data_loader = make_dataloader(dataset, batch_size=batch_size, shuffle=True)

        # minibatch is a list of (particles, tilt, indices)
        for i, minibatch in enumerate(data_loader):
            assert isinstance(minibatch, Sequence)
            assert len(minibatch) == 3

            # We have 100 particles. For all but the last iteration *
            # for all but the last iteration (100//7 = 14), we'll have 7 images each
            if i < (dataset.Np // batch_size):
                assert minibatch[0].shape == (batch_size * ntilts, dataset.D, dataset.D)
                assert minibatch[1].shape == (batch_size * ntilts,)
                assert minibatch[2].shape == (batch_size,)

            # and 100 % 7 = 2 for the very last one
            else:
                assert minibatch[0].shape == (
                    (dataset.Np % batch_size) * ntilts,
                    dataset.D,
                    dataset.D,
                )
                assert minibatch[1].shape == ((dataset.Np % batch_size) * ntilts,)
                assert minibatch[2].shape == (dataset.Np % batch_size,)


@pytest.mark.parametrize(
    "particles, indices",
    [
        ("hand", None),
        ("hand-tilt", None),
        ("toy.mrcs", None),
        pytest.param(
            "toy.mrcs",
            "random-100",
            marks=pytest.mark.xfail(raises=NotImplementedError),
        ),
        ("toy.txt", None),
        ("toy.star", None),
        ("tilts.star", None),
    ],
    indirect=True,
)
@pytest.mark.parametrize(
    "batch_size, buffer_size",
    [
        (5, 20),
        (10, 40),
        (27, 81),
        pytest.param(
            40, 10, marks=pytest.mark.xfail(reason="buffer must be 0 mod batch")
        ),
    ],
)
def test_data_shuffler(particles, indices, batch_size, buffer_size):
    ind = None if indices.path is None else load_pkl(indices.path)
    dataset = ImageDataset(mrcfile=particles.path, ind=ind)
    data_loader = DataShuffler(dataset, batch_size=batch_size, buffer_size=buffer_size)

    # minibatch is a list of (particles, tilt, indices)
    epoch1_indices, epoch2_indices = list(), list()
    for i, minibatch in enumerate(data_loader):
        assert isinstance(minibatch, Sequence)
        assert len(minibatch) == 3

        assert minibatch[0].shape == (batch_size, dataset.D, dataset.D)
        assert minibatch[1].shape == (batch_size,)
        assert minibatch[2].shape == (batch_size,)
        epoch1_indices.append(minibatch[2])

    for i, minibatch in enumerate(data_loader):
        epoch2_indices.append(minibatch[2])

    epoch1_indices = np.concatenate(epoch1_indices)
    epoch2_indices = np.concatenate(epoch2_indices)

    N = len(epoch1_indices)
    # epochs should have all the indices exactly once
    assert sorted(epoch1_indices) == list(range(N)), epoch1_indices
    assert sorted(epoch2_indices) == list(range(N)), epoch2_indices

    # epochs should be shuffled differently
    assert any(epoch1_indices != epoch2_indices), epoch1_indices
