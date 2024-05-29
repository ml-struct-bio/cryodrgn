import pytest
import os
import numpy as np
import torch
from cryodrgn import fft
from cryodrgn.source import ImageSource
from cryodrgn.lattice import Lattice


def test_shifted_image():
    torch.manual_seed(15321)
    imgs = ImageSource.from_file(os.path.join(pytest.DATADIR, "hand.mrcs")).images()
    img = imgs[0]
    D = img.shape[0]
    ht = fft.ht2_center(img)
    ht = fft.symmetrize_ht(ht)
    D += 1

    lattice = Lattice(D)
    ht = ht.view(1, -1)  # type: ignore
    trans = torch.tensor([5.0, 10.0]).view(1, 1, 2)
    ht_shifted = lattice.translate_ht(ht, trans)
    ht_np = ht_shifted.view(D, D)[0:-1, 0:-1]

    img_shifted = fft.ihtn_center(ht_np)
    new_arr = torch.Tensor(np.load(os.path.join(pytest.DATADIR, "im_shifted.npy")))
    assert torch.allclose(new_arr, img_shifted, atol=1e-4)
