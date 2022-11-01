import os.path

import numpy as np
import torch

from cryodrgn import fft, mrc
from cryodrgn.lattice import Lattice

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


def test_shifted_image():
    torch.manual_seed(15321)
    imgs, _ = mrc.parse_mrc(f"{DATA_FOLDER}/hand.mrcs")
    img = imgs[0]
    D = img.shape[0]
    ht = fft.ht2_center(img)
    ht = fft.symmetrize_ht(ht)
    D += 1

    lattice = Lattice(D)
    ht = torch.tensor(ht.astype(np.float32)).view(1, -1)

    trans = torch.tensor([5.0, 10.0]).view(1, 1, 2)
    ht_shifted = lattice.translate_ht(ht, trans)
    ht_np = ht_shifted.view(D, D).numpy()[0:-1, 0:-1]

    img_shifted = fft.ihtn_center(ht_np)
    assert np.allclose(np.load(f"{DATA_FOLDER}/im_shifted.npy"), img_shifted)
