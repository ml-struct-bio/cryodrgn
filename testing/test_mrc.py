# coding: utf-8
import sys, os
from cryodrgn import mrc
import numpy as np
data, _ = mrc.parse_mrc('{}/data/toy_projections.mrcs'.format(DIR), lazy=True)
data2, _ = mrc.parse_mrc('{}/data/toy_projections.mrcs'.format(DIR), lazy=False)
data1=np.asarray([x.get() for x in data])
assert (data1==data2).all()
print('ok')

from cryodrgn import dataset
data2 = dataset.load_particles('{}/data/toy_projections.star'.format(DIR))
assert (data1==data2).all()
print('ok')

data2 = dataset.load_particles('{}/data/toy_projections.txt'.format(DIR))
assert (data1==data2).all()
print('ok')

print('all ok')
