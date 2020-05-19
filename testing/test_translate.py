'''
'''
import numpy as np
import sys, os
import argparse
import pickle
import matplotlib.pyplot as plt 

import torch
import torch.nn as nn

sys.path.insert(0,'../lib-python')
import fft
import models
import mrc
from lattice import Lattice

imgs,_ = mrc.parse_mrc('data/hand.mrcs')
img = imgs[0]
D = img.shape[0]
ht = fft.ht2_center(img)
ht = fft.symmetrize_ht(ht)
D += 1

lattice = Lattice(D)
model = models.FTSliceDecoder(D**2, D, 10,10,nn.ReLU)

coords = lattice.coords[...,0:2]/2
ht = torch.tensor(ht.astype(np.float32)).view(1,-1)

trans = torch.tensor([5.,10.]).view(1,1,2)
ht_shifted = lattice.translate_ht(ht, trans)
ht_np = ht_shifted.view(D,D).numpy()[0:-1, 0:-1]

img_shifted = fft.ihtn_center(ht_np)

plt.figure()
plt.imshow(img)
plt.figure()
plt.imshow(img_shifted)
plt.show()
