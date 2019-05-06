'''
Train a NN to model a 3D EM density map
    given 2D projections with angular assignments
'''
import numpy as np
import sys, os
import argparse
import pickle
from datetime import datetime as dt

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader

sys.path.insert(0,os.path.abspath(os.path.dirname(__file__))+'/lib-python')
import mrc
import utils
import fft
import dataset

from lattice import Lattice
from models import FTSliceDecoder

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrcs)')
    parser.add_argument('rot', help='Rotation matrix for each particle (.pkl)')
    parser.add_argument('--trans', help='Optionally provide translations for each particle (.pkl)')
    parser.add_argument('-o', '--outdir', type=os.path.abspath, required=True, help='Output directory to save model')
    parser.add_argument('--load', type=os.path.abspath, help='Initialize training from a checkpoint')
    parser.add_argument('--checkpoint', type=int, default=1, help='Checkpointing interval in N_EPOCHS (default: %(default)s)')
    parser.add_argument('--log-interval', type=int, default=1000, help='Logging interval in N_IMGS (default: %(default)s)')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')
    parser.add_argument('--seed', type=int, default=np.random.randint(0,100000), help='Random seed')

    group = parser.add_argument_group('Training parameters')
    group.add_argument('-n', '--num-epochs', type=int, default=10, help='Number of training epochs (default: %(default)s)')
    group.add_argument('-b','--batch-size', type=int, default=100, help='Minibatch size (default: %(default)s)')
    group.add_argument('--wd', type=float, default=0, help='Weight decay in Adam optimizer (default: %(default)s)')
    group.add_argument('--lr', type=float, default=1e-3, help='Learning rate in Adam optimizer (default: %(default)s)')

    group = parser.add_argument_group('Network Architecture')
    group.add_argument('--layers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--dim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')

    return parser

def save_checkpoint(model, lattice, optim, epoch, norm, out_mrc, out_weights):
    model.eval()
    vol = model.eval_volume(lattice.coords, lattice.D, norm)
    mrc.write(out_mrc, vol.astype(np.float32))
    torch.save({
        'norm': norm,
        'epoch':epoch,
        'model_state_dict':model.state_dict(),
        'optimizer_state_dict':optim.state_dict(),
        }, out_weights)

def train(model, lattice, optim, y, rot, trans=None):
    model.train()
    B = y.size(0)
    D = lattice.D
    yt = model(lattice.coords @ rot)
    yt = yt.view(-1, D, D)
    if trans:
        y = model.translate_ht(lattice.coords[:,0:2]/2, y.view(B,-1), trans.unsqueeze(1))
        y = y.view(-1, D, D)
    loss = F.mse_loss(yt, y)
    loss.backward()
    optim.step()
    return loss.item()

def main(args):
    t1 = dt.now()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # set the random seed
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    ## set the device
    use_cuda = torch.cuda.is_available()
    device = torch.device('cuda' if use_cuda else 'cpu')
    log('Use cuda {}'.format(use_cuda))
    if use_cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)

    # load the particles
    data = dataset.MRCData(args.particles)
    D = data.D
    Nimg = data.N

    lattice = Lattice(D)
    model = FTSliceDecoder(3, D, args.layers, args.dim, nn.ReLU)

    optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

    if args.load:
        log('Loading model weights from {}'.format(args.load))
        checkpoint = torch.load(args.load)
        model.load_state_dict(checkpoint['model_state_dict'])
        #optim.load_state_dict(checkpoint['optimizer_state_dict'])

    rot = torch.tensor(utils.load_pkl(args.rot))
    if args.trans: trans = torch.tensor(utils.load_pkl(args.trans))

    data_generator = DataLoader(data, batch_size=args.batch_size, shuffle=True)
    for epoch in range(args.num_epochs):
        loss_accum = 0
        batch_it = 0
        for batch, ind in data_generator:
            batch_it += len(ind)
            y = batch.to(device)
            r = rot[ind]
            t = trans[ind] if args.trans else None
            loss_item = train(model, lattice, optim, batch.to(device), r, t)
            loss_accum += loss_item
            if batch_it % args.log_interval == 0:
                log('# [Train Epoch: {}/{}] [{}/{} images] loss={:.4f}'.format(epoch+1, args.num_epochs, batch_it, Nimg, loss_item))
        log('# =====> Epoch: {} Average loss = {:.4}'.format(epoch+1, loss_accum/Nimg))
        if args.checkpoint and epoch % args.checkpoint == 0:
            out_mrc = '{}/reconstruct.{}.mrc'.format(args.outdir,epoch)
            out_weights = '{}/weights.{}.pkl'.format(args.outdir,epoch)
            save_checkpoint(model, lattice, optim, epoch, data.norm, out_mrc, out_weights)

    ## save model weights and evaluate the model on 3D lattice
    out_mrc = '{}/reconstruct.mrc'.format(args.outdir)
    out_weights = '{}/weights.pkl'.format(args.outdir)
    save_checkpoint(model, lattice, optim, epoch, data.norm, out_mrc, out_weights)
   
    td = dt.now()-t1
    log('Finsihed in {} ({} per epoch)'.format(td, td/args.num_epochs))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)
