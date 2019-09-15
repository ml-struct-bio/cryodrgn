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
import lie_tools

from lattice import Lattice
from bnb import BNBHet
from models import HetOnlyVAE
from beta_schedule import get_beta_schedule, LinearSchedule
from losses import EquivarianceLoss

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrc)')
    parser.add_argument('-o', '--outdir', type=os.path.abspath, required=True, help='Output directory to save model')
    parser.add_argument('--load', type=os.path.abspath, help='Initialize training from a checkpoint')
    parser.add_argument('--checkpoint', type=int, default=1, help='Checkpointing interval in N_EPOCHS (default: %(default)s)')
    parser.add_argument('--log-interval', type=int, default=1000, help='Logging interval in N_IMGS (default: %(default)s)')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')
    parser.add_argument('--seed', type=int, default=np.random.randint(0,100000), help='Random seed')
    parser.add_argument('--invert-data', action='store_true', help='Invert data sign')

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
    group.add_argument('--norm', type=float, nargs=2, default=None, help='Data normalization as shift, 1/scale (default: mean, std of dataset)')

    group = parser.add_argument_group('Branch and bound')
    group.add_argument('--l-start', type=int,default=12, help='Starting L radius (default: %(default)s)')
    group.add_argument('--l-end', type=int, default=20, help='End L radius (default: %(default)s)')
    group.add_argument('--l-end-it',type=int,default=100000, help='default: %(default)s')
    group.add_argument('--rotate', action='store_true', help='Speedup BNB with image rotation')
    group.add_argument('--t-extent', type=float, default=5, help='+/- pixels to search over translations (default: %(default)s)')
    group.add_argument('--no-trans', action='store_true', help="Don't search over translations")
    group.add_argument('--bnb-start', type=int, default=1, help='Number of initial BNNB epochs (default: %(default)s)')

    group = parser.add_argument_group('Encoder Network')
    group.add_argument('--qlayers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--qdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--encode-mode', default='resid', choices=('conv','resid','mlp','tilt'), help='Type of encoder network (default: %(default)s)')
    group.add_argument('--zdim', type=int, default=1, help='Dimension of latent variable')
    group.add_argument('--enc-mask', type=int, help='Circulask mask of image for encoder')

    group = parser.add_argument_group('Decoder Network')
    group.add_argument('--players', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--pdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--enc-type', choices=('geom_ft','geom_full','geom_lowf','geom_nohighf','linear_lowf','none'), default='linear_lowf', help='Type of positional encoding')
    group.add_argument('--domain', choices=('hartley','fourier'), default='fourier')
    return parser

def pretrain(model, lattice, optim, minibatch, tilt):
    y, yt = minibatch
    use_tilt = yt is not None
    B = y.size(0)

    model.train()
    optim.zero_grad()

    rot = lie_tools.random_SO3(B, device=y.device)
    z = torch.randn((B,model.z_dim), device=y.device)

    # reconstruct circle of pixels instead of whole image
    mask = lattice.get_circular_mask(lattice.D//2)
    def gen_slice(R):
        return model.decode(lattice.coords[mask] @ R, z).view(B,-1)
    
    y = y.view(B,-1)[:, mask]
    if use_tilt:
        yt = yt.view(B,-1)[:, mask]
        gen_loss = .5*F.mse_loss(gen_slice(rot), y) + .5*F.mse_loss(gen_slice(tilt @ rot), yt)
    else:
        gen_loss = F.mse_loss(gen_slice(rot), y)
    
    gen_loss.backward()
    optim.step()
    return gen_loss.item()

def train(model, lattice, bnb, optim, minibatch, beta, beta_control=None, equivariance=None, rotated_images=None, enc_only=False, no_trans=False):
    y, yt = minibatch
    use_tilt = yt is not None
    D = lattice.D
    B = y.size(0)

    # VAE inference of z
    model.train()
    optim.zero_grad()
    input_ = (y,yt) if use_tilt else (y,)
    z_mu, z_logvar = model.encode(*input_)
    z = model.reparameterize(z_mu, z_logvar)

    if equivariance is not None:
        lamb, equivariance_loss = equivariance
        eq_loss = equivariance_loss(y, z_mu)

    # BNB inference of pose
    model.eval()
    with torch.no_grad():
        if no_trans:
            if rotated_images is None:
                rot = bnb.opt_theta(y, z, None if enc_only else yt)
            else:
                rot = bnb.opt_theta_rot(y, rotated_images, z)
        else:
            rot, trans = bnb.opt_theta_trans(y, z, None if enc_only else yt)

    model.train()

    # reconstruct circle of pixels instead of whole image
    mask = lattice.get_circular_mask(lattice.D//2)
    def gen_slice(R):
        return model.decode(lattice.coords[mask] @ R, z).view(B,-1)
    def translate(img):
        img = model.decoder.translate_ht(lattice.freqs2d[mask], img, trans.unsqueeze(1))
        return img.view(B,-1)

    y = y.view(B,-1)[:, mask]
    if use_tilt: yt = yt.view(B,-1)[:, mask]
    if not no_trans:
        y = translate(y)
        if use_tilt: yt = translate(yt)    

    if use_tilt:
        gen_loss = .5*F.mse_loss(gen_slice(rot), y) + .5*F.mse_loss(gen_slice(bnb.tilt @ rot), yt)
    else:
        gen_loss = F.mse_loss(gen_slice(rot), y)

    kld = -0.5 * torch.mean(1 + z_logvar - z_mu.pow(2) - z_logvar.exp())
    if torch.isnan(kld):
        log(z_mu[0])
        log(z_logvar[0])
        raise RuntimeError('KLD is nan')

    if beta_control is None:
        loss = gen_loss + beta*kld/mask.sum()
    else:
        loss = gen_loss + args.beta_control*(beta-kld)**2/mask.sum()

    if equivariance is not None:
        loss += lamb*eq_loss

    loss.backward()
    optim.step()
    save_pose = [rot.detach().cpu().numpy()]
    if not no_trans:
        save_pose.append(trans.detach().cpu().numpy())
    return gen_loss.item(), kld.item(), loss.item(), eq_loss.item() if equivariance else None, save_pose

def eval_z(model, lattice, data, batch_size, device, use_tilt=False):
    assert not model.training
    z_mu_all = []
    z_logvar_all = []
    data_generator = DataLoader(data, batch_size=batch_size, shuffle=False)
    for minibatch in data_generator:
        y = minibatch[0].to(device)
        if use_tilt: yt = minibatch[1].to(device)
        input_ = (y,yt) if use_tilt else (y,)
        z_mu, z_logvar = model.encode(*input_)
        z_mu_all.append(z_mu.detach().cpu().numpy())
        z_logvar_all.append(z_logvar.detach().cpu().numpy())
    z_mu_all = np.vstack(z_mu_all)
    z_logvar_all = np.vstack(z_logvar_all)
    return z_mu_all, z_logvar_all

def save_checkpoint(model, lattice, optim, epoch, norm, bnb_pose, z_mu, z_logvar, out_mrc_dir, out_weights, out_z):
    '''Save model weights, latent encoding z, and decoder volumes'''
    # save model weights
    torch.save({
        'epoch':epoch,
        'model_state_dict':model.state_dict(),
        'optimizer_state_dict':optim.state_dict(),
        'bnb_pose': bnb_pose
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
        data = dataset.MRCData(args.particles, keepreal=args.rotate, norm=args.norm, invert_data=args.invert_data)
        tilt = None
    else:
        assert args.encode_mode == 'tilt'
        data = dataset.TiltMRCData(args.particles, args.tilt, keepreal=args.rotate, norm=args.norm, invert_data=args.invert_data)
        tilt = torch.tensor(utils.xrot(args.tilt_deg).astype(np.float32))
    D = data.D
    Nimg = data.N

    lattice = Lattice(D, extent=0.5)
    if args.enc_mask: args.enc_mask = lattice.get_circular_mask(args.enc_mask)
    model = HetOnlyVAE(lattice, args.qlayers, args.qdim, args.players, args.pdim,
                args.zdim, encode_mode=args.encode_mode, enc_mask=args.enc_mask, enc_type=args.enc_type, domain=args.domain)
    bnb = BNBHet(model, lattice, args.l_start, args.l_end, tilt, args.t_extent)
    if args.rotate: 
        assert args.enc_only
        theta = torch.arange(1,12,dtype=torch.float32)*2*np.pi/12 # 11 angles 

    if args.equivariance:
        assert args.equivariance > 0, 'Regularization weight must be positive'
        equivariance_lambda = LinearSchedule(0, args.equivariance, 10000, args.equivariance_end_it)
        equivariance_loss = EquivarianceLoss(model, D)

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

    # training loop
    num_epochs = args.num_epochs
    data_iterator = DataLoader(data, batch_size=args.batch_size, shuffle=True)
    for epoch in range(start_epoch, num_epochs):
        if epoch < args.bnb_start:
            log('[Train Epoch: {}/{}] Using branch and no bound'.format(epoch+1, args.num_epochs))
        kld_accum = 0
        gen_loss_accum = 0
        loss_accum = 0
        eq_loss_accum = 0
        batch_it = 0 
        bnb_pose = []
        for batch in data_iterator:
            ind = batch[-1]
            batch = (batch[0].to(device), None) if tilt is None else (batch[0].to(device), batch[1].to(device))
            batch_it += len(batch[0])
            global_it = Nimg*epoch + batch_it

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
                loss = pretrain(model, lattice, optim, batch, bnb.tilt)
                gen_loss = kld = eq_loss = -1
            else:
                gen_loss, kld, loss, eq_loss, pose = train(model, lattice, bnb, optim, batch, beta, args.beta_control, equivariance_tuple, rotated_images=yr, enc_only=args.enc_only, no_trans=args.no_trans)
            
                # logging
                bnb_pose.append((ind.cpu().numpy(),pose))
                kld_accum += kld*len(ind)
                gen_loss_accum += gen_loss*len(ind)
                if args.equivariance:eq_loss_accum += eq_loss*len(ind)

            loss_accum += loss*len(ind)
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
                z_mu, z_logvar = eval_z(model, lattice, data, args.batch_size, device, tilt is not None)
                save_checkpoint(model, lattice, optim, epoch, data.norm, bnb_pose, z_mu, z_logvar, out_mrc, out_weights, out_z)

    ## save model weights and evaluate the model on 3D lattice
    model.eval()
    out_mrc = '{}/reconstruct'.format(args.outdir)
    out_weights = '{}/weights.pkl'.format(args.outdir)
    out_z = '{}/z.pkl'.format(args.outdir)
    with torch.no_grad():
        z_mu, z_logvar = eval_z(model, lattice, data, args.batch_size, device, tilt is not None)
        save_checkpoint(model, lattice, optim, epoch, data.norm, bnb_pose, z_mu, z_logvar, out_mrc, out_weights, out_z)
    
    td = dt.now()-t1
    log('Finsihed in {} ({} per epoch)'.format(td, td/(num_epochs-start_epoch)))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

