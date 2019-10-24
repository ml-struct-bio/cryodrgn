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
import lie_tools

from lattice import Lattice
from bnb import BNBHomo, BNBHomoRot
from models import FTPositionalDecoder, PositionalDecoder

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrcs)')
    parser.add_argument('-o', '--outdir', type=os.path.abspath, required=True, help='Output directory to save model')
    parser.add_argument('--norm', type=float, nargs=2, default=None, help='Data normalization as shift, 1/scale (default: mean, std of dataset)')
    parser.add_argument('--load', type=os.path.abspath, help='Initialize training from a checkpoint')
    parser.add_argument('--checkpoint', type=int, default=1, help='Checkpointing interval in N_EPOCHS (default: %(default)s)')
    parser.add_argument('--log-interval', type=int, default=1000, help='Logging interval in N_IMGS (default: %(default)s)')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')
    parser.add_argument('--seed', type=int, default=np.random.randint(0,100000), help='Random seed')
    parser.add_argument('--invert-data', action='store_true', help='Invert data sign')
    parser.add_argument('--window', action='store_true', help='Real space windowing of dataset')

    group = parser.add_argument_group('Tilt series')
    group.add_argument('--tilt', help='Particle stack file (.mrcs)')
    group.add_argument('--tilt-deg', type=float, default=45, help='X-axis tilt offset in degrees (default: %(default)s)')

    group = parser.add_argument_group('Training parameters')
    group.add_argument('--t-extent', type=float, default=5, help='+/- pixels to search over translations')
    group.add_argument('--t-ngrid', type=float, default=7, help='Initial grid size for translations')
    group.add_argument('--no-trans', action='store_true', help="Don't search over translations")
    group.add_argument('--pretrain', type=int, default=10000, help='Number of initial iterations with random poses (default: %(default)s)')
    group.add_argument('--bnb-freq', type=int, default=1, help='Frequency of pose inference (default: every %(default)s epochs)')
    group.add_argument('-n', '--num-epochs', type=int, default=10, help='Number of training epochs (default: %(default)s)')
    group.add_argument('-b','--batch-size', type=int, default=10, help='Minibatch size (default: %(default)s)')
    group.add_argument('--wd', type=float, default=0, help='Weight decay in Adam optimizer (default: %(default)s)')
    group.add_argument('--lr', type=float, default=1e-4, help='Learning rate in Adam optimizer (default: %(default)s)')
    group.add_argument('--l-start', type=int,default=12, help='Starting L radius (default: %(default)s)')
    group.add_argument('--l-end', type=int, default=20, help='End L radius (default: %(default)s)')
    group.add_argument('--l-end-it',type=int,default=100000, help='default: %(default)s')
    group.add_argument('--probabilistic', action='store_true', help='Use probabilistic bound')

    group = parser.add_argument_group('Network Architecture')
    group.add_argument('--layers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--dim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--pe-type', choices=('geom_ft','geom_full','geom_lowf','geom_nohighf','linear_lowf'), default='geom_lowf', help='Type of positional encoding')
    group.add_argument('--domain', choices=('hartley','fourier'), default='hartley')

    return parser

def save_checkpoint(model, lattice, bnb_pose, optim, epoch, norm, out_mrc, out_weights, out_poses):
    model.eval()
    vol = model.eval_volume(lattice.coords, lattice.D, lattice.extent, norm)
    mrc.write(out_mrc, vol.astype(np.float32))
    torch.save({
        'norm': norm,
        'epoch':epoch,
        'model_state_dict':model.state_dict(),
        'optimizer_state_dict':optim.state_dict(),
        }, out_weights)
    with open(out_poses,'wb') as f:
        pickle.dump(bnb_pose, f)

def pretrain(model, lattice, optim, batch, tilt=None):
    y, yt = batch
    B = y.size(0)
    model.train()
    optim.zero_grad()

    mask = lattice.get_circular_mask(lattice.D//2)
    def gen_slice(R):
        slice_ = model(lattice.coords[mask] @ R)
        return slice_.view(B,-1)

    rot = lie_tools.random_SO3(B, device=y.device)

    y = y.view(B,-1)[:, mask]
    if tilt is not None:
        yt = yt.view(B,-1)[:, mask]
        loss = .5*F.mse_loss(gen_slice(rot), y) + .5*F.mse_loss(gen_slice(tilt @ rot), yt)
    else:
        loss = F.mse_loss(gen_slice(rot), y) 
    loss.backward()
    optim.step()
    return loss.item()

def sort_bnb_poses(bnb_pose):
    ind = [x[0] for x in bnb_pose]
    ind = np.concatenate(ind)
    rot = [x[1][0] for x in bnb_pose]
    rot = np.concatenate(rot)
    rot = rot[np.argsort(ind)]
    if len(bnb_pose[0][1]) == 2:
        trans = [x[1][1] for x in bnb_pose]
        trans = np.concatenate(trans)
        trans = trans[np.argsort(ind)]
        return (rot,trans)
    return (rot,)

def train(model, lattice, bnb, optim, batch, tilt=None, no_trans=False, poses=None):
    y, yt = batch
    B = y.size(0)

    # pose inference
    if poses is not None:
        rot = poses[0]
        if not no_trans: trans = poses[1]
    else: # BNB
        model.eval()
        with torch.no_grad():
            if no_trans:
                rot = bnb.opt_theta(y, yt)
            else:
                rot, trans = bnb.opt_theta_trans(y, yt)

    # reconstruct circle of pixels instead of whole image
    mask = lattice.get_circular_mask(lattice.D//2)
    def gen_slice(R):
        slice_ = model(lattice.coords[mask] @ R)
        return slice_.view(B,-1)
    def translate(img):
        img = lattice.translate_ht(img, trans.unsqueeze(1), mask)
        return img.view(B,-1)

    # Train model 
    model.train()
    optim.zero_grad()

    y = y.view(B,-1)[:, mask]
    if tilt is not None: yt = yt.view(B,-1)[:, mask]
    if not no_trans:
        y = translate(y)
        if tilt is not None: yt = translate(yt)

    if tilt is not None:
        loss = .5*F.mse_loss(gen_slice(rot), y) + .5*F.mse_loss(gen_slice(tilt @ rot), yt)
    else:
        loss = F.mse_loss(gen_slice(rot), y) 
    loss.backward()
    optim.step()
    save_pose = [rot.detach().cpu().numpy()]
    if not no_trans:
        save_pose.append(trans.detach().cpu().numpy())
    return loss.item(), save_pose

def main(args):
    log(args)
    if args.probabilistic:
        assert args.no_trans, 'Probabilistic bound not implemented with translations yet'
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
        data = dataset.MRCData(args.particles, norm=args.norm, invert_data=args.invert_data, window=args.window)
        tilt = None
    else:
        data = dataset.TiltMRCData(args.particles, args.tilt, norm=args.norm, invert_data=args.invert_data, window=args.window)
        tilt = torch.tensor(utils.xrot(args.tilt_deg).astype(np.float32))
    D = data.D
    Nimg = data.N

    lattice = Lattice(D, extent=0.5)
    if args.domain == 'fourier':
        model = FTPositionalDecoder(3, D, args.layers, args.dim, nn.ReLU, enc_type=args.pe_type)
    else:
        model = PositionalDecoder(3, D, args.layers, args.dim, nn.ReLU, enc_type=args.pe_type)

    if args.no_trans:
        bnb = BNBHomoRot(model, lattice, args.l_start, args.l_end, tilt, args.probabilistic)
    else:    
        bnb = BNBHomo(model, lattice, args.l_start, args.l_end, tilt, t_extent=args.t_extent, t_ngrid=args.t_ngrid)
    log(model)
    log('{} parameters in model'.format(sum(p.numel() for p in model.parameters() if p.requires_grad)))
    optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

    if args.load:
        args.pretrain = 0
        log('Loading checkpoint from {}'.format(args.load))
        checkpoint = torch.load(args.load)
        model.load_state_dict(checkpoint['model_state_dict'])
        optim.load_state_dict(checkpoint['optimizer_state_dict'])
        start_epoch = checkpoint['epoch']+1
        assert args.num_epochs > start_epoch
        model.train()
    else:
        start_epoch = 0

    data_iterator = DataLoader(data, batch_size=args.batch_size, shuffle=True)

    # pretrain decoder with random poses
    global_it = 0
    pretrain_epoch = 0
    log('Using random poses for {} iterations'.format(args.pretrain))
    while global_it < args.pretrain:
        for batch in data_iterator:
            global_it += len(batch[0])
            batch = (batch[0].to(device), None) if tilt is None else (batch[0].to(device), batch[1].to(device))
            loss = pretrain(model, lattice, optim, batch, tilt=bnb.tilt)
            if global_it % args.log_interval == 0:
                log(f'[Pretrain Iteration {global_it}] loss={loss:4f}')
            if global_it > args.pretrain:
                break

    # training loop
    for epoch in range(start_epoch, args.num_epochs):
        t2 = dt.now()
        batch_it = 0
        loss_accum = 0
        bnb_pose = []
        if epoch % args.bnb_freq != 0:
            log('Using previous iteration poses')
        for batch in data_iterator:
            ind = batch[-1]
            batch = (batch[0].to(device), None) if tilt is None else (batch[0].to(device), batch[1].to(device))
            batch_it += len(batch[0])
            
            # train the model
            if epoch % args.bnb_freq != 0:
                p = [torch.tensor(x[ind]) for x in sorted_poses]
            else: 
                p = None
            loss_item, pose = train(model, lattice, bnb, optim, batch, tilt, args.no_trans, poses=p) 
            bnb_pose.append((ind.cpu().numpy(),pose))
           
            # logging
            loss_accum += loss_item*len(batch[0])
            if batch_it % args.log_interval == 0:
                log('# [Train Epoch: {}/{}] [{}/{} images] loss={:.4f}'.format(epoch+1, args.num_epochs, batch_it, Nimg, loss_item))
        log('# =====> Epoch: {} Average loss = {:.4}; Finished in {}'.format(epoch+1, loss_accum/Nimg, dt.now() - t2))

        # sort bnb_pose
        sorted_poses = sort_bnb_poses(bnb_pose) if bnb_pose else None

        if args.checkpoint and epoch % args.checkpoint == 0:
            out_mrc = '{}/reconstruct.{}.mrc'.format(args.outdir,epoch)
            out_weights = '{}/weights.{}.pkl'.format(args.outdir,epoch)
            out_poses = '{}/pose.{}.pkl'.format(args.outdir, epoch)
            save_checkpoint(model, lattice, sorted_poses, optim, epoch, data.norm, out_mrc, out_weights, out_poses)

    ## save model weights and evaluate the model on 3D lattice
    out_mrc = '{}/reconstruct.mrc'.format(args.outdir)
    out_weights = '{}/weights.pkl'.format(args.outdir)
    out_poses = '{}/pose.pkl'.format(args.outdir, epoch)
    save_checkpoint(model, lattice, sorted_poses, optim, epoch, data.norm, out_mrc, out_weights, out_poses)
   
    td = dt.now()-t1
    log('Finsihed in {} ({} per epoch)'.format(td, td/(args.num_epochs-start_epoch)))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

