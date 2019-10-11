'''
Train a VAE for heterogeneous reconstruction with known pose
'''
import numpy as np
import sys, os
import argparse
import pickle
from datetime import datetime as dt
import json

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader

sys.path.insert(0,os.path.abspath(os.path.dirname(__file__))+'/lib-python')
import mrc
import utils
import fft
import lie_tools
import dataset
import ctf

from models import HetOnlyVAE
from lattice import Lattice
from beta_schedule import get_beta_schedule, LinearSchedule
from losses import EquivarianceLoss

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('particles', type=os.path.abspath, help='Particles (.mrcs)')
    parser.add_argument('-o', '--outdir', type=os.path.abspath, required=True, help='Output directory to save model')
    parser.add_argument('--ctf', metavar='pkl', type=os.path.abspath, help='CTF parameters (.pkl)')
    parser.add_argument('--poses', type=os.path.abspath, nargs='*', required=True, help='Image rotations and optionally translations (.pkl)')
    parser.add_argument('--tscale', type=float, default=1.0, help='Scale translations by this amount')
    parser.add_argument('--load', type=os.path.abspath, help='Initialize training from a checkpoint')
    parser.add_argument('--checkpoint', type=int, default=1, help='Checkpointing interval in N_EPOCHS (default: %(default)s)')
    parser.add_argument('--log-interval', type=int, default=1000, help='Logging interval in N_IMGS (default: %(default)s)')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')
    parser.add_argument('--seed', type=int, default=np.random.randint(0,100000), help='Random seed')
    parser.add_argument('--invert-data', action='store_true', help='Invert data sign')
    parser.add_argument('--no-window', dest='window', action='store_false', help='Do not window dataset')
    parser.add_argument('--ind', type=os.path.abspath, help='Filter indices')

    group = parser.add_argument_group('Tilt series')
    group.add_argument('--tilt', help='Particles (.mrcs)')
    group.add_argument('--tilt-deg', type=float, default=45, help='X-axis tilt offset in degrees (default: %(default)s)')

    group = parser.add_argument_group('Training parameters')
    group.add_argument('-n', '--num-epochs', type=int, default=10, help='Number of training epochs (default: %(default)s)')
    group.add_argument('-b','--batch-size', type=int, default=100, help='Minibatch size (default: %(default)s)')
    group.add_argument('--wd', type=float, default=0, help='Weight decay in Adam optimizer (default: %(default)s)')
    group.add_argument('--lr', type=float, default=1e-4, help='Learning rate in Adam optimizer (default: %(default)s)')
    group.add_argument('--beta', default=1.0, help='Choice of beta schedule or a constant for KLD weight (default: %(default)s)')
    group.add_argument('--beta-control', type=float, help='KL-Controlled VAE gamma. Beta is KL target. (default: %(default)s)')
    group.add_argument('--equivariance', type=float, help='Strength of equivariance loss (default: %(default)s)')
    group.add_argument('--equivariance-end-it', type=int, default=100000, help='It at which equivariance max (default: %(default)s)')
    group.add_argument('--norm', type=float, nargs=2, default=None, help='Data normalization as shift, 1/scale (default: mean, std of dataset)')

    group = parser.add_argument_group('Encoder Network')
    group.add_argument('--qlayers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--qdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--encode-mode', default='resid', choices=('conv','resid','mlp','tilt'), help='Type of encoder network (default: %(default)s)')
    group.add_argument('--zdim', type=int, default=1, help='Dimension of latent variable')
    group.add_argument('--enc-mask', type=int, help='Circular mask of image for encoder')
    group.add_argument('--use-real', action='store_true', help='Use real space image for encoder')

    group = parser.add_argument_group('Decoder Network')
    group.add_argument('--players', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--pdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--pe-type', choices=('geom_ft','geom_full','geom_lowf','geom_nohighf','linear_lowf','none'), default='linear_lowf', help='Type of positional encoding')
    group.add_argument('--domain', choices=('hartley','fourier'), default='fourier')
    return parser

def train(model, lattice, y, yt, rot, trans, optim, beta, beta_control=None, equivariance=None, tilt=None, ctf_params=None, yr=None):
    use_tilt = yt is not None
    use_ctf = ctf_params is not None
    model.train()
    optim.zero_grad()
    B = y.size(0)
    D = lattice.D
    if use_ctf:
        freqs = lattice.freqs2d.unsqueeze(0).expand(B,*lattice.freqs2d.shape)/ctf_params[:,0].view(B,1,1)
        c = ctf.compute_ctf(freqs, *torch.split(ctf_params[:,1:], 1, 1)).view(B,D,D)
    
    # translate the image
    if trans is not None:
        y = lattice.translate_ht(y.view(B,-1), trans.unsqueeze(1)).view(B,D,D)
        if yt is not None: yt = lattice.translate_ht(yt.view(B,-1), trans.unsqueeze(1)).view(B,D,D)

    # inference of z
    if yr is not None:
        input_ = (yr,)
    else:
        input_ = (y,yt) if yt is not None else (y,)
        if use_ctf: input_ = (x*c.sign() for x in input_) # phase flip by the ctf
    z_mu, z_logvar = model.encode(*input_)
    z = model.reparameterize(z_mu, z_logvar)

    # decode 
    mask = lattice.get_circular_mask(D//2) # restrict to circular mask
    y_recon = model.decode(lattice.coords[mask]/lattice.extent/2 @ rot, z).view(B,-1)
    if use_ctf: y_recon *= c.view(B,-1)[:,mask]
    gen_loss = F.mse_loss(y_recon, y.view(B,-1)[:, mask])

    # decode the tilt series
    if use_tilt:
        y_recon_tilt = model.decode(lattice.coords[mask]/lattice.extent/2 @ tilt @ rot, z)
        if use_ctf: y_recon_tilt *= c.view(B,-1)[:,mask]
        gen_loss = .5*gen_loss + .5*F.mse_loss(y_recon_tilt, yt.view(B,-1)[:,mask])

    # latent loss
    kld = -0.5 * torch.mean(1 + z_logvar - z_mu.pow(2) - z_logvar.exp())

    if beta_control is None:
        loss = gen_loss + beta*kld/mask.sum()
    else:
        loss = gen_loss + args.beta_control*(beta-kld)**2/mask.sum()

    # extra equivariance loss term
    if equivariance is not None:
        lamb, equivariance_loss = equivariance
        eq_loss = equivariance_loss(y, z_mu)
        loss += lamb*eq_loss

    loss.backward()
    optim.step()
    return gen_loss.item(), kld.item(), loss.item(), eq_loss.item() if equivariance else None

def eval_z(model, lattice, data, batch_size, device, trans=None, use_tilt=False, ctf_params=None, use_real=False):
    assert not model.training
    z_mu_all = []
    z_logvar_all = []
    data_generator = DataLoader(data, batch_size=batch_size, shuffle=False)
    for minibatch in data_generator:
        ind = minibatch[-1]
        y = minibatch[0].to(device)
        yt = minibatch[1].to(device) if use_tilt else None
        B = len(ind)
        D = lattice.D
        if ctf_params is not None:
            freqs = lattice.freqs2d.unsqueeze(0).expand(B,*lattice.freqs2d.shape)/ctf_params[ind,0].view(B,1,1)
            c = ctf.compute_ctf(freqs, *torch.split(ctf_params[ind,1:], 1, 1)).view(B,D,D)
        if trans is not None:
            y = lattice.translate_ht(y.view(B,-1), trans[ind].unsqueeze(1)).view(B,D,D)
            if yt is not None: yt = lattice.translate_ht(yt.view(B,-1), trans[ind].unsqueeze(1)).view(B,D,D)
        if use_real:
            input_ = (torch.from_numpy(data.particles_real[ind]).to(device),)
        else:
            input_ = (y,yt) if yt is not None else (y,)
        if ctf_params is not None: 
            assert not use_real, "Not implemented"
            input_ = (x*c.sign() for x in input_) # phase flip by the ctf
        z_mu, z_logvar = model.encode(*input_)
        z_mu_all.append(z_mu.detach().cpu().numpy())
        z_logvar_all.append(z_logvar.detach().cpu().numpy())
    z_mu_all = np.vstack(z_mu_all)
    z_logvar_all = np.vstack(z_logvar_all)
    return z_mu_all, z_logvar_all
    
def save_checkpoint(model, lattice, optim, epoch, norm, z_mu, z_logvar, out_mrc_dir, out_weights, out_z):
    '''Save model weights, latent encoding z, and decoder volumes'''
    # save model weights
    torch.save({
        'epoch':epoch,
        'model_state_dict':model.state_dict(),
        'optimizer_state_dict':optim.state_dict(),
        }, out_weights)
    # save z
    with open(out_z,'wb') as f:
        pickle.dump(z_mu, f)
        pickle.dump(z_logvar, f)
    # save single structure at mean of z_mu
    vol = model.decoder.eval_volume(lattice.coords, lattice.D, lattice.extent, norm, z_mu.mean(axis=0))
    mrc.write(out_mrc_dir+'.mrc', vol.astype(np.float32))
    log('Saved {} with z = {}'.format(out_mrc_dir+'.mrc', z_mu.mean(axis=0)))
    # save trajectory of structures if zdim = 1
    if z_mu.shape[1] == 1:
        if not os.path.exists(out_mrc_dir): os.mkdir(out_mrc_dir)
        for i in range(5):
            pct = 10+i*20
            zz = np.percentile(z_mu, pct, keepdims=True)
            vol = model.decoder.eval_volume(lattice.coords, lattice.D, lattice.extent, norm, zz)
            mrc.write('{}/traj{}.mrc'.format(out_mrc_dir,int(pct)),vol)
            log('Saved {}/traj{}.mrc with z = {}'.format(out_mrc_dir, int(pct), zz))

def save_config(args, dataset, lattice, model, out_config):
    dataset_args = dict(particles=args.particles,
                        norm=dataset.norm,
                        invert_data=args.invert_data,
                        ind=args.ind,
                        keepreal=args.use_real,
                        window=args.window)
    if args.tilt is not None:
        dataset_args['particles_tilt'] = args.tilt
    lattice_args = dict(D=lattice.D,
                        extent=lattice.extent,
                        ignore_DC=lattice.ignore_DC)
    model_args = dict(qlayers=args.qlayers,
                      qdim=args.qdim,
                      players=args.players,
                      pdim=args.pdim,
                      zdim=args.zdim,
                      encode_mode=args.encode_mode,
                      enc_mask=args.enc_mask,
                      pe_type=args.pe_type,
                      domain=args.domain)
    with open(out_config,'wb') as f:
        pickle.dump(dict(dataset_args=dataset_args,
                         lattice_args=lattice_args,
                         model_args=model_args), f)

def main(args):
    log(args)
    t1 = dt.now()
    if args.outdir is not None and not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # set the random seed
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    # set the device
    use_cuda = torch.cuda.is_available()
    device = torch.device('cuda' if use_cuda else 'cpu')
    log('Use cuda {}'.format(use_cuda))
    if use_cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)

    # set beta schedule
    try:
        args.beta = float(args.beta)
    except ValueError: 
        assert args.beta_control, "Need to set beta control weight for schedule {}".format(args.beta)
    beta_schedule = get_beta_schedule(args.beta)

    # load the particles
    if args.ind is not None: 
        log('Filtering image dataset with {}'.format(args.ind))
        ind = pickle.load(open(args.ind,'rb'))
    else: ind = None
    if args.tilt is None:
        if args.encode_mode == 'conv':
            args.use_real = True
        data = dataset.MRCData(args.particles, norm=args.norm, invert_data=args.invert_data, ind=ind, keepreal=args.use_real, window=args.window)
        tilt = None
    else:
        assert args.encode_mode == 'tilt'
        data = dataset.TiltMRCData(args.particles, args.tilt, norm=args.norm, invert_data=args.invert_data, ind=ind, window=args.window, keepreal=args.use_real)
        tilt = torch.tensor(utils.xrot(args.tilt_deg).astype(np.float32))
    Nimg = data.N
    D = data.D

    if args.encode_mode == 'conv':
        assert D-1 == 64, "Image size must be 64x64 for convolutional encoder"

    # load poses
    assert len(args.poses) in (1,2)
    if len(args.poses) == 2: # rotation pickle, translation pickle
        poses = (utils.load_pkl(args.poses[0]), utils.load_pkl(args.poses[1]))
    else: # rotation pickle or poses pickle
        poses = utils.load_pkl(args.poses[0])
        if type(poses) != tuple: poses = (poses,)
    rots = torch.tensor(poses[0]).float()
    if args.ind is not None: rots = rots[ind]
    assert rots.shape == (Nimg,3,3)
    if len(args.poses) == 2:
        trans = args.tscale * torch.tensor(poses[1]).float()
        if args.ind is not None: trans = trans[ind]
        assert trans.shape == (Nimg,2)
    else: trans = None

    # load ctf
    if args.ctf is not None:
        if args.use_real:
            raise NotImplementedError("Not implemented with real-space encoder. Use phase-flipped images instead")
        log('Loading ctf params from {}'.format(args.ctf))
        ctf_params = utils.load_pkl(args.ctf)
        if args.ind is not None: ctf_params = ctf_params[ind]
        assert ctf_params.shape == (Nimg, 7)
        ctf.print_ctf_params(ctf_params[0])
        ctf_params = torch.tensor(ctf_params)
    else: ctf_params = None

    # instantiate model
    lattice = Lattice(D, extent=0.5)
    if args.enc_mask: 
        enc_mask = lattice.get_circular_mask(args.enc_mask)
        in_dim = enc_mask.sum()
    else:
        enc_mask = None
        in_dim = lattice.D**2 if not args.use_real else (lattice.D-1)**2
    model = HetOnlyVAE(lattice, args.qlayers, args.qdim, args.players, args.pdim,
                in_dim, args.zdim, encode_mode=args.encode_mode, enc_mask=enc_mask,
                enc_type=args.pe_type, domain=args.domain)
    log(model)
    log('{} parameters in model'.format(sum(p.numel() for p in model.parameters() if p.requires_grad)))

    # save configuration
    out_config = '{}/config.pkl'.format(args.outdir)
    save_config(args, data, lattice, model, out_config)

    if args.equivariance:
        assert args.equivariance > 0, 'Regularization weight must be positive'
        equivariance_lambda = LinearSchedule(0, args.equivariance, 10000, args.equivariance_end_it)
        equivariance_loss = EquivarianceLoss(model, D)

    optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

    # restart from checkpoint
    if args.load:
        log('Loading checkpoint from {}'.format(args.load))
        checkpoint = torch.load(args.load)
        model.load_state_dict(checkpoint['model_state_dict'])
        optim.load_state_dict(checkpoint['optimizer_state_dict'])
        start_epoch = checkpoint['epoch']+1
        model.train()
    else:
        start_epoch = 0

    # training loop
    data_generator = DataLoader(data, batch_size=args.batch_size, shuffle=True)
    num_epochs = args.num_epochs
    for epoch in range(start_epoch, num_epochs):
        gen_loss_accum = 0
        loss_accum = 0
        kld_accum = 0
        eq_loss_accum = 0
        batch_it = 0 
        for minibatch in data_generator:
            ind = minibatch[-1]
            y = minibatch[0].to(device)
            yt = minibatch[1].to(device) if tilt is not None else None
            B = len(ind)
            batch_it += B
            global_it = Nimg*epoch + batch_it

            beta = beta_schedule(global_it)
            if args.equivariance:
                lamb = equivariance_lambda(global_it)
                equivariance_tuple = (lamb, equivariance_loss)
            else:
                equivariance_tuple = None
           
            yr = torch.from_numpy(data.particles_real[ind]).to(device) if args.use_real else None
            rot = rots[ind]
            tran = trans[ind] if trans is not None else None
            ctf_param = ctf_params[ind] if ctf_params is not None else None
            gen_loss, kld, loss, eq_loss = train(model, lattice, y, yt, rot, tran, optim, beta, args.beta_control, equivariance_tuple, tilt, ctf_params=ctf_param, yr=yr)

            # logging
            gen_loss_accum += gen_loss*B
            kld_accum += kld*B
            loss_accum += loss*B
            if args.equivariance: eq_loss_accum += eq_loss*B

            if batch_it % args.log_interval == 0:
                eq_log = 'equivariance={:.4f}, lambda={:.4f}, '.format(eq_loss, lamb) if args.equivariance else ''
                log('# [Train Epoch: {}/{}] [{}/{} images] gen loss={:.4f}, kld={:.4f}, beta={:.4f}, {}loss={:.4f}'.format(epoch+1, num_epochs, batch_it, Nimg, gen_loss, kld, beta, eq_log, loss))
        eq_log = 'equivariance = {:.4f}, '.format(eq_loss_accum/Nimg) if args.equivariance else ''
        log('# =====> Epoch: {} Average gen loss = {:.4}, KLD = {:.4f}, {}total loss = {:.4f}'.format(epoch+1, gen_loss_accum/Nimg, kld_accum/Nimg, eq_log, loss_accum/Nimg))

        if args.checkpoint and epoch % args.checkpoint == 0:
            out_mrc = '{}/reconstruct.{}'.format(args.outdir,epoch)
            out_weights = '{}/weights.{}.pkl'.format(args.outdir,epoch)
            out_z = '{}/z.{}.pkl'.format(args.outdir, epoch)
            model.eval()
            with torch.no_grad():
                z_mu, z_logvar = eval_z(model, lattice, data, args.batch_size, device, trans, tilt is not None, ctf_params, args.use_real)
                save_checkpoint(model, lattice, optim, epoch, data.norm, z_mu, z_logvar, out_mrc, out_weights, out_z)

    # save model weights, latent encoding, and evaluate the model on 3D lattice
    out_mrc = '{}/reconstruct'.format(args.outdir)
    out_weights = '{}/weights.pkl'.format(args.outdir)
    out_z = '{}/z.pkl'.format(args.outdir)
    model.eval()
    with torch.no_grad():
        z_mu, z_logvar = eval_z(model, lattice, data, args.batch_size, device, trans, tilt is not None, ctf_params, args.use_real)
        save_checkpoint(model, lattice, optim, epoch, data.norm, z_mu, z_logvar, out_mrc, out_weights, out_z)
    
    td = dt.now()-t1
    log('Finsihed in {} ({} per epoch)'.format(td, td/(num_epochs-start_epoch)))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

