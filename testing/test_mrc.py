# coding: utf-8
import sys, os
DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,'{}/../lib-python'.format(DIR))
import mrc
import numpy as np
data, _, _ = mrc.parse_mrc('{}/data/toy_projections.mrcs'.format(DIR), lazy=True)
data2, _, _ = mrc.parse_mrc('{}/data/toy_projections.mrcs'.format(DIR), lazy=False)
data1=np.asarray([x.get() for x in data])
assert (data1==data2).all()
print('all ok')
