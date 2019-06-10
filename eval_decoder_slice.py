'''
Evaluate slices out of the decoder model
'''
import numpy as np
import sys, os
import argparse
import pickle
from datetime import datetime as dt
import matplotlib.pyplot as plt 

import torch
import torch.nn as nn

sys.path.insert(0,os.path.abspath(os.path.dirname(__file__))+'/lib-python')
import mrc
import utils
import fft
from lattice import Lattice
from models import FTSliceDecoder

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('weights', help='Model weights')
    parser.add_argument('rot', help='Orientations (.pkl)')
    parser.add_argument('-D', type=int, help='Image dimension D')
    parser.add_argument('--norm', type=float, nargs=2, default=None, help='Data normalization as shift, 1/scale (default: mean, std of dataset)')
    parser.add_argument('--first', type=int, help='First N rotations')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output mrcs')
    parser.add_argument('-v','--verbose',action='store_true',help='Increase verbosity')

    group = parser.add_argument_group('Architecture parameters')
    group.add_argument('-b','--batch-size', type=int, default=100, help='Minibatch size (default: %(default)s)')
    group.add_argument('--layers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--dim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    return parser

def main(args):
    log(args)
    t1 = dt.now()

    ## set the device
    use_cuda = torch.cuda.is_available()
    device = torch.device('cuda' if use_cuda else 'cpu')
    log('Use cuda {}'.format(use_cuda))
    if use_cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)

    assert args.D % 2 == 0
    D = args.D + 1
    lattice = Lattice(D)
    model = FTSliceDecoder(3,D,args.layers,args.dim,nn.ReLU)

    log('Loading weights from {}'.format(args.weights))
    checkpoint = torch.load(args.weights)
    model.load_state_dict(checkpoint['model_state_dict'])

    rots = utils.load_pkl(args.rot)
    if args.first: rots = rots[:args.first]
    log('Evaluating {} slices'.format(len(rots)))

    imgs = []
    model.eval()
    for rot in np.array_split(rots,len(rots)//args.batch_size + 1):
        rot = torch.tensor(rot)
        y_recon = model(lattice.coords @ rot)
        y_recon = y_recon.view(-1,D,D)
        imgs.append(y_recon.detach().cpu().numpy())
        break

    imgs = np.vstack(imgs)
    imgs = imgs[:,0:-1,0:-1]
    imgs = np.asarray([fft.ihtn_center(x) for x in imgs])
    log(imgs.shape)
    mrc.write(args.o, imgs)

    td = dt.now()-t1
    log('Finsihed in {}'.format(td))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

