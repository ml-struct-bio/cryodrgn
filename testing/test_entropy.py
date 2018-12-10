import numpy as np
import sys, os
import argparse
import pickle
from datetime import datetime as dt

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
from torch.distributions import Normal

#sys.path.insert(0,os.path.abspath(os.path.dirname(__file__))+'/lib-python')
#import mrc
#import utils
#import fft

sys.path.insert(0,'..')
import vae_rot

avg = []
std = torch.tensor([2.3407,1.0999,1.2962])
for _ in range(100):
    w_eps = torch.randn_like(std)*std
    e = vae_rot.so3_entropy(w_eps, std)
    avg.append(e)
    print(e)

w_eps = torch.tensor([-.46,-1.54,-1.96])
e = vae_rot.so3_entropy(w_eps, std)
print(e)
avg.append(e)

a = torch.mean(torch.Tensor(avg))
print('average: {}'.format(a))
print(np.log(8*np.pi**2))
