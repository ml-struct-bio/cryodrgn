import pytest
import os
import numpy as np
from typing import Sequence
from torch.utils.data.sampler import BatchSampler, RandomSampler
from torch.utils.data import DataLoader
from cryodrgn.dataset import DataShuffler, ImageDataset, make_dataloader


@pytest.mark.parametrize("particles", ["hand"], indirect=True)
def test_particles(particles):
    # compare data for some random indices against data generated from known
    # correct values after normalization/fft/symmetrization etc.
    dataset = ImageDataset(mrcfile=particles, invert_data=True)
    indices = np.array([43, 12, 53, 64, 31, 56, 75, 63, 27, 62, 96])
    particles_arr, _, _ = dataset[indices]

    # the tolerance (atol) of 1e-5 accounts for minor fft differences between
    # numpy fft (legacy) and torch fft
    assert np.allclose(
        np.load(os.path.join(pytest.data_dir, "hand_11_particles.npy")),
        particles_arr.cpu().numpy(),
        atol=1e-5,
    )


@pytest.mark.parametrize(
    "particles",
    ["hand", "hand-tilt", "toy.mrcs", "toy.txt", "toy.star", "tilts.star"],
    indirect=True,
)
class TestLoading:
    @pytest.mark.parametrize("batch_size", [7, 20, 43])
    def test_loading_slow(self, particles, batch_size):
        dataset = ImageDataset(mrcfile=particles)
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

    @pytest.mark.parametrize("batch_size", [25, 61])
    def test_loading_fast(self, particles, batch_size):
        dataset = ImageDataset(mrcfile=particles)

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

    def test_data_shuffler(self, particles):
        dataset = ImageDataset(mrcfile=particles)
        data_loader = DataShuffler(dataset, batch_size=5, buffer_size=20)
        epoch1_indices, epoch2_indices = [], []

        # minibatch is a list of (particles, tilt, indices)
        for i, minibatch in enumerate(data_loader):
            assert isinstance(minibatch, Sequence)
            assert len(minibatch) == 3

            assert minibatch[0].shape == (5, dataset.D, dataset.D)
            assert minibatch[1].shape == (5,)
            assert minibatch[2].shape == (5,)
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
