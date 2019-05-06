'''
Evaluate the encoder trained by bnb_het.py
'''
import numpy as np
import sys, os
import argparse
import pickle
from datetime import datetime as dt
import matplotlib.pyplot as plt 

import torch
from torch.utils.data import DataLoader

sys.path.insert(0,os.path.abspath(os.path.dirname(__file__))+'/lib-python')
import mrc
import utils
import dataset
from lattice import Lattice
from models import HetOnlyVAE

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrcs)')
    parser.add_argument('weights', help='Model weights')
    parser.add_argument('--tilt', help='Particle stack file for tilt series (.mrcs')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output pickle')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')

    group = parser.add_argument_group('Architecture parameters')
    group.add_argument('-b','--batch-size', type=int, default=100, help='Minibatch size (default: %(default)s)')
    group.add_argument('--qlayers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--qdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--zdim', type=int, default=1, help='Dimension of latent variable')
    group.add_argument('--encode-mode', default='resid', choices=('conv','resid','mlp','tilt'), help='Type of encoder network (default: %(default)s)')
    group.add_argument('--players', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--pdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
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

    # load the particles
    if args.tilt:
        data = dataset.TiltMRCData(args.particles,args.tilt)
    else:
        data = dataset.MRCData(args.particles)

    Nimg, D = data.N, data.D
    lattice = Lattice(D)
    model = HetOnlyVAE(lattice, D*D, args.qlayers, args.qdim, args.players, args.pdim,
                args.zdim, encode_mode=args.encode_mode)

    log('Loading weights from {}'.format(args.weights))
    checkpoint = torch.load(args.weights)
    model.load_state_dict(checkpoint['model_state_dict'])

    model.eval()
    z_all = []
    data_iterator = DataLoader(data, batch_size=args.batch_size, shuffle=False)
    for batch in data_iterator:
        B = len(batch[0])
        y = batch[0].to(device)
        if args.tilt:
            yt = batch[0].to(device)
            mu, logvar = model.encode((y.view(B,-1), yt.view(B,-1)))
        else:
            mu, logvar = model.encode(y.view(B,-1))
        z_all.append(mu.detach().cpu().numpy())

    if args.o:
        z_all = np.vstack(z_all)
        log(np.mean(z_all))
        with open(args.o,'wb') as f:
            pickle.dump(z_all, f)

    td = dt.now()-t1
    log('Finsihed in {}'.format(td))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

