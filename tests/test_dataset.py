import os.path
from torch.utils.data.sampler import BatchSampler, RandomSampler
from torch.utils.data import DataLoader
from cryodrgn.dataset import ImageDataset

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


def test_loading_slow():
    dataset = ImageDataset(mrcfile=f"{DATA_FOLDER}/toy_projections.mrcs")
    # We could use the following, but this will result in 7 independent calls to
    # the underlying ImageDataset's __getitem__ - not efficient
    data_loader = DataLoader(dataset, batch_size=7, shuffle=True)

    for i, minibatch in enumerate(data_loader):
        assert isinstance(minibatch, list)
        assert len(minibatch) == 3  # minibatch is a list of (particles, tilt, indices)

        # We have 1000 particles. For all but the last iteration *
        # for all but the last iteration (1000//7 = 142), we'll have 7 images each
        if i < 142:
            assert minibatch[0].shape == (7, 1, 31, 31)
            assert minibatch[1].shape == (7, 1, 31, 31)
            assert minibatch[2].shape == (7,)
        # and 1000 % 7 = 6 for the very last one
        else:
            assert minibatch[0].shape == (6, 1, 31, 31)
            assert minibatch[1].shape == (6, 1, 31, 31)
            assert minibatch[2].shape == (6,)


def test_loading_fast():
    dataset = ImageDataset(mrcfile=f"{DATA_FOLDER}/toy_projections.mrcs")
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

        # We have 1000 particles. For all but the last iteration *
        # for all but the last iteration (1000//7 = 142), we'll have 7 images each
        if i < 142:
            assert minibatch[0].shape == (7, 31, 31)
            assert minibatch[1].shape == (7, 31, 31)
            assert minibatch[2].shape == (7,)
        # and 1000 % 7 = 6 for the very last one
        else:
            assert minibatch[0].shape == (6, 31, 31)
            assert minibatch[1].shape == (6, 31, 31)
            assert minibatch[2].shape == (6,)
