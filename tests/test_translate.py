import os.path
import numpy as np
import torch
from cryodrgn import fft
from cryodrgn.source import ImageSource
from cryodrgn.lattice import Lattice

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


def test_shifted_image():
    torch.manual_seed(15321)
    imgs = ImageSource.from_file(f"{DATA_FOLDER}/hand.mrcs").images()
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
    assert torch.allclose(
        torch.Tensor(np.load(f"{DATA_FOLDER}/im_shifted.npy")), img_shifted, atol=1e-4
    )
