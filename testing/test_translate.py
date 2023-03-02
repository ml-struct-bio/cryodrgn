import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import cryodrgn.fft
import cryodrgn.models
import cryodrgn.mrc
from cryodrgn.lattice import Lattice
from cryodrgn.source import ImageSource

imgs = ImageSource.from_file("data/hand.mrcs").images()
img = imgs[0]
D = img.shape[0]
ht = cryodrgn.fft.ht2_center(img)
ht = cryodrgn.fft.symmetrize_ht(ht)
D += 1

lattice = Lattice(D)
model = cryodrgn.models.FTSliceDecoder(D**2, D, 10, 10, nn.ReLU)

coords = lattice.coords[..., 0:2] / 2
ht = torch.tensor(ht).view(1, -1)

trans = torch.tensor([5.0, 10.0]).view(1, 1, 2)
ht_shifted = lattice.translate_ht(ht, trans)
ht_np = ht_shifted.view(D, D)[0:-1, 0:-1]

img_shifted = cryodrgn.fft.ihtn_center(ht_np)

plt.figure()
plt.imshow(img)
plt.figure()
plt.imshow(img_shifted)
plt.show()
