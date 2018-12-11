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

sys.path.insert(0,'../lib-python')
import lie_tools

avg = []
std = torch.tensor([2.3407,1.0999,1.2962])
for _ in range(10):
    w_eps = torch.randn_like(std)*std
    e = lie_tools.so3_entropy_old(w_eps, std)
    avg.append(e)
    print(e)

w_eps = torch.tensor([-.46,-1.54,-1.96])
e = lie_tools.so3_entropy_old(w_eps, std)
print(e)
avg.append(e)

# test new multi sample so3_entropy
w_eps = torch.tensor([-.46,-1.54,-1.96])
w_eps2 = torch.tensor([1.,1.,1.])
std = torch.tensor([2.3407,1.0999,1.2962])
std2 = torch.tensor([1.,1.,1.])
e = lie_tools.so3_entropy_old(w_eps2, std2)
print(e)
e = lie_tools.so3_entropy(torch.stack([w_eps,w_eps2]), torch.stack([std,std2]))
print(e)

a = torch.mean(torch.Tensor(avg))
print('average: {}'.format(a))
print(np.log(8*np.pi**2))
