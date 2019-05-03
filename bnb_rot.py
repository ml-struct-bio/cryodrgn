'''
Homogeneous NN reconstruction with BNB optimization of orientation
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
from bnb import BNNBHomo, BNBHomo, BNBHomoRot
from models import FTSliceDecoder
from losses import EquivarianceLoss
from beta_schedule import LinearSchedule

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrcs)')
    parser.add_argument('-o', '--outdir', type=os.path.abspath, required=True, help='Output directory to save model')
    parser.add_argument('--load', type=os.path.abspath, help='Initialize training from a checkpoint')
    parser.add_argument('--checkpoint', type=int, default=1, help='Checkpointing interval in N_EPOCHS (default: %(default)s)')
    parser.add_argument('--log-interval', type=int, default=1000, help='Logging interval in N_IMGS (default: %(default)s)')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')
    parser.add_argument('--seed', type=int, default=np.random.randint(0,100000), help='Random seed')

    group = parser.add_argument_group('Tilt series')
    group.add_argument('--tilt', help='Particle stack file (.mrcs)')
    group.add_argument('--tilt-deg', type=float, default=45, help='X-axis tilt offset in degrees (default: %(default)s)')

    group = parser.add_argument_group('Training parameters')
    group.add_argument('--t-extent', type=float, default=5, help='+/- pixels to search over translations')
    group.add_argument('--no-trans', action='store_true', help="Don't search over translations")
    group.add_argument('-n', '--num-epochs', type=int, default=10, help='Number of training epochs (default: %(default)s)')
    group.add_argument('-b','--batch-size', type=int, default=100, help='Minibatch size (default: %(default)s)')
    group.add_argument('--wd', type=float, default=0, help='Weight decay in Adam optimizer (default: %(default)s)')
    group.add_argument('--lr', type=float, default=1e-3, help='Learning rate in Adam optimizer (default: %(default)s)')
    group.add_argument('--l-start', type=int,default=12, help='Starting L radius (default: %(default)s)')
    group.add_argument('--l-end', type=int, default=20, help='End L radius (default: %(default)s)')
    group.add_argument('--l-end-it',type=int,default=100000, help='default: %(default)s')

    group = parser.add_argument_group('Network Architecture')
    group.add_argument('--layers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--dim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')

    return parser

def save_checkpoint(model, lattice, bnb_pose, optim, epoch, norm, out_mrc, out_weights):
    model.eval()
    vol = model.eval_volume(lattice.coords, lattice.D, norm)
    mrc.write(out_mrc, vol.astype(np.float32))
    torch.save({
        'norm': norm,
        'epoch':epoch,
        'model_state_dict':model.state_dict(),
        'optimizer_state_dict':optim.state_dict(),
        'bnb_pose': bnb_pose
        }, out_weights)

def train(model, lattice, bnb, optim, batch, L, tilt=None, no_trans=False):
    y, yt = batch
    B = y.size(0)
    # BNB inference of orientation 
    model.eval()
    with torch.no_grad():
        if no_trans:
            rot = bnb.opt_theta(y, L, yt)
        else:
            rot, trans = bnb.opt_theta_trans(y, L, yt)
            vlog(trans[0])
    # Train model 
    model.train()
    optim.zero_grad()
    y_recon = model(lattice.coords @ rot)
    y_recon = y_recon.view(-1, lattice.D, lattice.D)
    if not no_trans:
        y = model.translate_ht(lattice.coords[:,0:2]/2, y.view(B,-1), trans.unsqueeze(1))
        y = y.view(-1, lattice.D, lattice.D)
    loss = F.mse_loss(y_recon,y)
    if tilt is not None:
        y_recon_tilt = model(lattice.coords @ tilt @ rot)
        y_recon_tilt = y_recon_tilt.view(-1, lattice.D, lattice.D)
        if not no_trans:
            yt = model.translate_ht(lattice.coords[:,0:2]/2, yt.view(B,-1), trans.unsqueeze(1))
            yt = yt.view(-1, lattice.D, lattice.D)
        loss = .5*loss + .5*F.mse_loss(y_recon_tilt,yt)
    loss.backward()
    optim.step()
    save_pose = [rot.detach().cpu().numpy()]
    if not no_trans:
        save_pose.append(trans.detach().cpu().numpy())
    return loss.item(), save_pose

def main(args):
    log(args)
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
    if args.tilt is None:
        data = dataset.MRCData(args.particles)
        tilt = None
    else:
        data = dataset.TiltMRCData(args.particles, args.tilt)
        theta = args.tilt_deg*np.pi/180
        tilt = np.array([[1.,0.,0.],
                        [0, np.cos(theta), -np.sin(theta)],
                        [0, np.sin(theta), np.cos(theta)]]).astype(np.float32)
        tilt = torch.tensor(tilt)
    D = data.D
    Nimg = data.N

    lattice = Lattice(D)
    model = FTSliceDecoder(3, D, args.layers, args.dim, nn.ReLU)
    bnnb = BNNBHomo(model, lattice, tilt, t_extent=args.t_extent)
    if args.no_trans:
        bnb = BNBHomoRot(model, lattice, args.l_start, args.l_end, tilt)
    else:    
        bnb = BNBHomo(model, lattice, args.l_start, args.l_end, tilt, t_extent=args.t_extent)

    optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

    if args.load:
        log('Loading checkpoint from {}'.format(args.load))
        checkpoint = torch.load(args.load)
        model.load_state_dict(checkpoint['model_state_dict'])
        optim.load_state_dict(checkpoint['optimizer_state_dict'])
        start_epoch = checkpoint['epoch']+1
        assert args.num_epochs > start_epoch
        model.train()
    else:
        start_epoch = 0

    if args.l_start == -1:
        Lsched = lambda x: None
    else:
        Lsched = LinearSchedule(args.l_start,args.l_end,0,args.l_end_it)

    # training loop
    data_generator = DataLoader(data, batch_size=args.batch_size, shuffle=True)
    for epoch in range(start_epoch, args.num_epochs):
        batch_it = 0
        loss_accum = 0
        bnb_pose = []
        for batch in data_generator:
            ind = batch[-1]
            batch = (batch[0].to(device), None) if tilt is None else (batch[0].to(device), batch[1].to(device))
            batch_it += len(batch[0])
            global_it = Nimg*epoch+batch_it
            L = Lsched(global_it)
            if L: L = int(L)
            
            # train the model
            if epoch < 1:
                loss_item, pose = train(model, lattice, bnnb, optim, batch, L, tilt, args.no_trans)
            else:
                L = None
                loss_item, pose = train(model, lattice, bnb, optim, batch, L, tilt, args.no_trans)
           
            # logging
            bnb_pose.append((ind.cpu().numpy(),pose))
            loss_accum += loss_item*len(batch[0])
            if batch_it % args.log_interval == 0:
                log('# [Train Epoch: {}/{}] [{}/{} images] loss={:.4f}'.format(epoch+1, args.num_epochs, batch_it, Nimg, loss_item))
        log('# =====> Epoch: {} Average loss = {:.4}'.format(epoch+1, loss_accum/Nimg))

        if args.checkpoint and epoch % args.checkpoint == 0:
            out_mrc = '{}/reconstruct.{}.mrc'.format(args.outdir,epoch)
            out_weights = '{}/weights.{}.pkl'.format(args.outdir,epoch)
            save_checkpoint(model, lattice, bnb_pose, optim, epoch, data.norm, out_mrc, out_weights)

    ## save model weights and evaluate the model on 3D lattice
    out_mrc = '{}/reconstruct.mrc'.format(args.outdir)
    out_weights = '{}/weights.pkl'.format(args.outdir)
    save_checkpoint(model, lattice, bnb_pose, optim, epoch, data.norm, out_mrc, out_weights)
   
    td = dt.now()-t1
    log('Finsihed in {} ({} per epoch)'.format(td, td/args.num_epochs))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

