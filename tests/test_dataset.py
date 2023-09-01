import os.path
import numpy as np
from torch.utils.data.sampler import BatchSampler, RandomSampler
from torch.utils.data import DataLoader
from cryodrgn.dataset import DataShuffler, ImageDataset

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


def test_particles():
    dataset = ImageDataset(mrcfile=f"{DATA_FOLDER}/hand.mrcs", invert_data=True)
    # Compare data for some random indices against data generated from known correct values
    # after normalization/fft/symmetrization etc.
    indices = np.array([43, 12, 53, 64, 31, 56, 75, 63, 27, 62, 96])
    particles, _, _ = dataset[indices]
    # The atol of 1e-5 accounts for minor fft differences between numpy fft (legacy) and torch fft
    assert np.allclose(
        np.load(f"{DATA_FOLDER}/hand_11_particles.npy"),
        particles.cpu().numpy(),
        atol=1e-5,
    )


def test_loading_slow():
    dataset = ImageDataset(mrcfile=f"{DATA_FOLDER}/hand.mrcs")
    # We could use the following, but this will result in 7 independent calls to
    # the underlying ImageDataset's __getitem__ - not efficient
    data_loader = DataLoader(dataset, batch_size=7, shuffle=True)

    for i, minibatch in enumerate(data_loader):
        assert isinstance(minibatch, list)
        assert len(minibatch) == 3  # minibatch is a list of (particles, tilt, indices)

        # We have 100 particles. For all but the last iteration *
        # for all but the last iteration (100//7 = 14), we'll have 7 images each
        if i < 14:
            assert minibatch[0].shape == (7, 65, 65)
            assert minibatch[1].shape == (7, 65, 65)
            assert minibatch[2].shape == (7,)
        # and 100 % 7 = 2 for the very last one
        else:
            assert minibatch[0].shape == (2, 65, 65)
            assert minibatch[1].shape == (2, 65, 65)
            assert minibatch[2].shape == (2,)


def test_loading_fast():
    dataset = ImageDataset(mrcfile=f"{DATA_FOLDER}/hand.mrcs")
    # A faster way to load is to use BatchSampler with RandomSampler
    # see https://stackoverflow.com/questions/61458305
    data_loader = DataLoader(
        dataset,
        sampler=BatchSampler(
            RandomSampler(dataset),
            batch_size=7,
            drop_last=False,  # Or SequentialSampler
        ),
        batch_size=None,  # Do not let torch add a new batch dimension at the front
    )

    for i, minibatch in enumerate(data_loader):
        assert isinstance(minibatch, list)
        assert len(minibatch) == 3  # minibatch is a list of (particles, tilt, indices)

        # We have 100 particles. For all but the last iteration *
        # for all but the last iteration (100//7 = 14), we'll have 7 images each
        if i < 14:
            assert minibatch[0].shape == (7, 65, 65)
            assert minibatch[1].shape == (7, 65, 65)
            assert minibatch[2].shape == (7,)
        # and 100 % 7 = 2 for the very last one
        else:
            assert minibatch[0].shape == (2, 65, 65)
            assert minibatch[1].shape == (2, 65, 65)
            assert minibatch[2].shape == (2,)


def test_data_shuffler():
    dataset = ImageDataset(mrcfile=f"{DATA_FOLDER}/hand.mrcs")
    # We could use the following, but this will result in 7 independent calls to
    # the underlying ImageDataset's __getitem__ - not efficient
    data_loader = DataShuffler(dataset, batch_size=5, buffer_size=20)

    epoch1_indices, epoch2_indices = [], []

    for i, minibatch in enumerate(data_loader):
        assert len(minibatch) == 3  # minibatch is a list of (particles, tilt, indices)

        # We have 100 particles. For all but the last iteration *
        # for all but the last iteration (100//7 = 14), we'll have 7 images each
        assert minibatch[0].shape == (5, 65, 65)
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
    assert any(epoch1_indices != epoch2_indices), epoch1_indices  # Should be reshuffled
