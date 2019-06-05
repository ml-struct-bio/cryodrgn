'''
Heterogeneous NN reconstruction with BNB optimization of orientation
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

from lattice import EvenLattice
from bnb import BNNBHet, BNBHet
from models import HetOnlyVAE
from beta_schedule import get_beta_schedule, LinearSchedule
from losses import EquivarianceLoss

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrc)')
    parser.add_argument('-o', '--outdir', type=os.path.abspath, required=True, help='Output directory to save model')
    parser.add_argument('--norm', type=float, nargs=2, default=None, help='Data normalization as shift, 1/scale (default: mean, std of dataset)')
    parser.add_argument('--load', type=os.path.abspath, help='Initialize training from a checkpoint')
    parser.add_argument('--checkpoint', type=int, default=1, help='Checkpointing interval in N_EPOCHS (default: %(default)s)')
    parser.add_argument('--log-interval', type=int, default=1000, help='Logging interval in N_IMGS (default: %(default)s)')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')
    parser.add_argument('--seed', type=int, default=np.random.randint(0,100000), help='Random seed')

    group = parser.add_argument_group('Tilt series')
    group.add_argument('--tilt', help='Particle stack file (.mrcs)')
    group.add_argument('--tilt-deg', type=float, default=45, help='X-axis tilt offset in degrees (default: %(default)s)')
    group.add_argument('--enc-only', action='store_true', help='Use the tilt pair only in VAE and not in BNB search')
    
    group = parser.add_argument_group('Training parameters')
    group.add_argument('-n', '--num-epochs', type=int, default=10, help='Number of training epochs (default: %(default)s)')
    group.add_argument('-b','--batch-size', type=int, default=100, help='Minibatch size (default: %(default)s)')
    group.add_argument('--wd', type=float, default=0, help='Weight decay in Adam optimizer (default: %(default)s)')
    group.add_argument('--lr', type=float, default=1e-3, help='Learning rate in Adam optimizer (default: %(default)s)')
    group.add_argument('--beta', default=1.0, help='Choice of beta schedule or a constant for KLD weight (default: %(default)s)')
    group.add_argument('--beta-control', type=float, help='KL-Controlled VAE gamma. Beta is KL target. (default: %(default)s)')
    group.add_argument('--equivariance', type=float, help='Strength of equivariance loss (default: %(default)s)')
    group.add_argument('--equivariance-end-it', type=int, default=100000, help='It at which equivariance max (default: %(default)s)')

    group = parser.add_argument_group('Branch and bound')
    group.add_argument('--l-start', type=int,default=12, help='Starting L radius (default: %(default)s)')
    group.add_argument('--l-end', type=int, default=20, help='End L radius (default: %(default)s)')
    group.add_argument('--l-end-it',type=int,default=100000, help='default: %(default)s')
    group.add_argument('--rotate', action='store_true', help='Speedup BNB with image rotation')
    group.add_argument('--t-extent', type=float, default=5, help='+/- pixels to search over translations')
    group.add_argument('--no-trans', action='store_true', help="Don't search over translations")
    group.add_argument('--bnb-start', type=int, default=1, help='Number of initial BNNB epochs')

    group = parser.add_argument_group('Encoder Network')
    group.add_argument('--qlayers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--qdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--encode-mode', default='resid', choices=('conv','resid','mlp','tilt'), help='Type of encoder network (default: %(default)s)')
    group.add_argument('--zdim', type=int, default=1, help='Dimension of latent variable')
    group.add_argument('--enc-mask', type=int, help='Circulask mask of image for encoder')

    group = parser.add_argument_group('Decoder Network')
    group.add_argument('--players', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--pdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    return parser

def eval_volume(model, lattice, D, zval, norm):
    '''Evaluate the model on a nz x ny x nx lattice'''
    zdim = len(zval)
    z = torch.zeros(D**2,zdim, dtype=torch.float32)
    z += torch.tensor(zval, dtype=torch.float32)

    vol_f = np.zeros((D,D,D),dtype=np.float32)
    assert not model.training
    # evaluate the volume by zslice to avoid memory overflows
    for i, dz in enumerate(np.linspace(-1,1,D,endpoint=True)):
        x = lattice.coords + torch.tensor([0,0,dz], dtype=torch.float32)
        x = torch.cat((x,z),dim=-1)
        with torch.no_grad():
            y = model.decoder.decode(x)
            y = y[...,0] - y[...,1]
            y = y.view(D, D).cpu().numpy()
        vol_f[i] = y
    vol_f = vol_f*norm[1] + norm[0]
    vol = fft.ihtn_center(vol_f[0:-1,0:-1,0:-1])
    return vol, vol_f

def train(model, lattice, bnb, optim, minibatch, L, beta, beta_control=None, equivariance=None, rotated_images=None, enc_only=False, no_trans=False):
    y, yt = minibatch
    B = y.size(0)
    model.train()
    optim.zero_grad()

    # inference of z
    D = lattice.D
    input_ = (y,) if yt is None else (y,yt)
    z_mu, z_logvar = model.encode(*input_)
    z = model.reparameterize(z_mu, z_logvar)

    # inference of pose
    model.eval()
    if no_trans:
        if rotated_images is None:
            rot = bnb.opt_theta(y, z, None if enc_only else yt, L=L)
        else:
            rot = bnb.opt_theta_rot(y, rotated_images, z, L=L)
    else:
        rot, trans = bnb.opt_theta_trans(y, z, None if enc_only else yt, L=L)
    model.train()

    y_recon = model.decode(rot, z)
    y_recon = y_recon.view(-1, D, D)
    if not no_trans:
        y = model.decoder.translate_ht(lattice.coords[:,0:2]/2, y.view(B,-1), trans.unsqueeze(1))
        y = y.view(-1, D, D)
    gen_loss = F.mse_loss(y_recon, y)
    if yt is not None: 
        y_recon_tilt = model.decode(bnb.tilt @ rot, z)
        y_recon_tilt = y_recon_tilt.view(-1, D, D)
        if not no_trans:
            yt = model.decoder.translate_ht(lattice.coords[:,0:2]/2, yt.view(B,-1), trans.unsqueeze(1))
            yt = yt.view(-1, lattice.D, lattice.D)
        gen_loss = .5*gen_loss + .5*F.mse_loss(y_recon_tilt, yt)

    kld = -0.5 * torch.mean(1 + z_logvar - z_mu.pow(2) - z_logvar.exp())

    if beta_control is None:
        loss = gen_loss + beta*kld/(D*D)
    else:
        loss = gen_loss + args.beta_control*(beta-kld)**2/(D*D)

    if equivariance is not None:
        lamb, equivariance_loss = equivariance
        eq_loss = equivariance_loss(y, z_mu)
        loss += lamb*eq_loss
    loss.backward()
    optim.step()
    save_pose = [rot.detach().cpu().numpy()]
    if not no_trans:
        save_pose.append(trans.detach().cpu().numpy())
    return gen_loss.item(), kld.item(), loss.item(), eq_loss.item() if equivariance else None, z_mu.detach().cpu().numpy().mean(0), save_pose

def save_checkpoint(model, lattice, z, bnb_pose, optim, epoch, norm, out_mrc, out_weights):
    model.eval()
    vol, vol_f = eval_volume(model, lattice, lattice.D, z, norm)
    mrc.write(out_mrc, vol.astype(np.float32))
    torch.save({
        'norm': norm,
        'epoch':epoch,
        'model_state_dict':model.state_dict(),
        'optimizer_state_dict':optim.state_dict(),
        'bnb_pose': bnb_pose
        }, out_weights)

def main(args):
    assert not args.rotate, "Not implemented with new BNB"

    log(args)
    t1 = dt.now()
    if args.outdir is not None and not os.path.exists(args.outdir):
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

    ## set beta schedule
    try:
        args.beta = float(args.beta)
    except ValueError: 
        assert args.beta_control, "Need to set beta control weight for schedule {}".format(args.beta)
    beta_schedule = get_beta_schedule(args.beta)

    # load the particles
    if args.tilt is None:
        data = dataset.MRCData(args.particles, keepreal=args.rotate, norm=args.norm)
        tilt = None
    else:
        assert args.encode_mode == 'tilt'
        data = dataset.TiltMRCData(args.particles, args.tilt, keepreal=args.rotate, norm=args.norm)
        theta = args.tilt_deg*np.pi/180
        tilt = np.array([[1.,0.,0.],
                        [0, np.cos(theta), -np.sin(theta)],
                        [0, np.sin(theta), np.cos(theta)]]).astype(np.float32)
        tilt = torch.tensor(tilt)
    D = data.D
    Nimg = data.N

    lattice = EvenLattice(D)
    if args.enc_mask: args.enc_mask = lattice.get_circular_mask(args.enc_mask)
    model = HetOnlyVAE(lattice, args.qlayers, args.qdim, args.players, args.pdim,
                args.zdim, encode_mode=args.encode_mode, enc_mask=args.enc_mask)
    bnnb = BNNBHet(model, lattice, tilt)
    bnb = BNBHet(model, lattice, args.l_start, args.l_end, tilt, args.t_extent)
    if args.rotate: 
        assert args.enc_only
        theta = torch.arange(1,12,dtype=torch.float32)*2*np.pi/12 # 11 angles 

    if args.equivariance:
        assert args.equivariance > 0, 'Regularization weight must be positive'
        equivariance_lambda = LinearSchedule(0, args.equivariance, 10000, args.equivariance_end_it)
        equivariance_loss = EquivarianceLoss(model,D,D)

    optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

    if args.load:
        log('Loading checkpoint from {}'.format(args.load))
        checkpoint = torch.load(args.load)
        model.load_state_dict(checkpoint['model_state_dict'])
        optim.load_state_dict(checkpoint['optimizer_state_dict'])
        start_epoch = checkpoint['epoch']+1
        model.train()
    else:
        start_epoch = 0

    if args.l_start == -1:
        Lsched = lambda x: None
    else:
        Lsched = LinearSchedule(args.l_start,args.l_end,0,args.l_end_it)

    # training loop
    num_epochs = args.num_epochs
    data_iterator = DataLoader(data, batch_size=args.batch_size, shuffle=True)
    for epoch in range(start_epoch, num_epochs):
        kld_accum = 0
        gen_loss_accum = 0
        loss_accum = 0
        eq_loss_accum = 0
        batch_it = 0 
        z_accum = np.zeros(args.zdim)
        bnb_pose = []
        for batch in data_iterator:
            ind = batch[-1]
            batch = (batch[0].to(device), None) if tilt is None else (batch[0].to(device), batch[1].to(device))
            batch_it += len(batch[0])
            global_it = Nimg*epoch + batch_it

            L = Lsched(global_it)
            if L: L = int(L)
            beta = beta_schedule(global_it)
            if args.equivariance:
                lamb = equivariance_lambda(global_it)
                equivariance_tuple = (lamb, equivariance_loss)
            else: equivariance_tuple = None

            if args.rotate:
                yr = torch.from_numpy(data.particles_real[ind]).to(device)
                yr = lattice.rotate(yr, theta)
                yr = fft.ht2_center(yr)
                yr = (yr-data.norm[0])/data.norm[1]
                yr = torch.from_numpy(yr.astype(np.float32)).to(device)
            else: yr = None

            # train the model
            if epoch < args.bnb_start:
                gen_loss, kld, loss, eq_loss, z, pose = train(model, lattice, bnnb, optim, batch, L, beta, args.beta_control, equivariance_tuple, rotated_images=yr, enc_only=args.enc_only, no_trans=True)
            else:
                L = None
                gen_loss, kld, loss, eq_loss, z, pose = train(model, lattice, bnb, optim, batch, L, beta, args.beta_control, equivariance_tuple, rotated_images=yr, enc_only=args.enc_only, no_trans=args.no_trans)
            
            # logging
            bnb_pose.append((ind.cpu().numpy(),pose))
            z_accum += z*len(ind)
            kld_accum += kld*len(ind)
            gen_loss_accum += gen_loss*len(ind)
            loss_accum += loss*len(ind)
            if args.equivariance:eq_loss_accum += eq_loss*len(ind)

            if batch_it % args.log_interval == 0:
                eq_log = 'equivariance={:.4f}, lambda={:.4f}, '.format(eq_loss, lamb) if args.equivariance else ''
                log('# [Train Epoch: {}/{}] [{}/{} images] gen loss={:.4f}, kld={:.4f}, beta={:.4f}, {}loss={:.4f}'.format(epoch+1, num_epochs, batch_it, Nimg, gen_loss, kld, beta, eq_log, loss))

        eq_log = 'equivariance = {:.4f}, '.format(eq_loss_accum/Nimg) if args.equivariance else ''
        log('# =====> Epoch: {} Average gen loss = {:.4}, KLD = {:.4f}, {}total loss = {:.4f}'.format(epoch+1, gen_loss_accum/Nimg, kld_accum/Nimg, eq_log, loss_accum/Nimg))

        if args.checkpoint and epoch % args.checkpoint == 0:
            out_mrc = '{}/reconstruct.{}.mrc'.format(args.outdir,epoch)
            out_weights = '{}/weights.{}.pkl'.format(args.outdir,epoch)
            save_checkpoint(model, lattice, z_accum/Nimg, bnb_pose, optim, epoch, data.norm, out_mrc, out_weights)

    ## save model weights and evaluate the model on 3D lattice
    model.eval()
    out_mrc = '{}/reconstruct.mrc'.format(args.outdir)
    out_weights = '{}/weights.pkl'.format(args.outdir)
    save_checkpoint(model, lattice, z_accum/Nimg, bnb_pose, optim, epoch, data.norm, out_mrc, out_weights)
    
    td = dt.now()-t1
    log('Finsihed in {} ({} per epoch)'.format(td, td/num_epochs))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

