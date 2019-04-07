'''
NN reconstruction with direct optimization of orientation
'''
import numpy as np
import sys, os
import argparse
import pickle
from datetime import datetime as dt

import torch
#torch.backends.cudnn.enabled=False
#torch.backends.cudnn.benchmark = True
import torch.nn as nn
import torch.nn.functional as F

sys.path.insert(0,os.path.abspath(os.path.dirname(__file__))+'/lib-python')
import mrc
import utils
import fft
import lie_tools
from models import BNBHetOpt, HetVAE, Lattice
from beta_schedule import get_beta_schedule, LinearSchedule
from losses import EquivarianceLoss

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrc)')
    parser.add_argument('-o', '--outdir', type=os.path.abspath, required=True, help='Output directory to save model')
    parser.add_argument('--load', type=os.path.abspath, help='Initialize training from a checkpoint')
    parser.add_argument('--checkpoint', type=int, default=5, help='Checkpointing interval in N_EPOCHS (default: %(default)s)')
    parser.add_argument('--log-interval', type=int, default=1000, help='Logging interval in N_IMGS (default: %(default)s)')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')
    parser.add_argument('--seed', type=int, default=np.random.randint(0,100000), help='Random seed')

    group = parser.add_argument_group('Tilt series')
    group.add_argument('--tilt', help='Particle stack file (.mrcs)')
    group.add_argument('--tilt-deg', type=float, default=45, help='X-axis tilt offset in degrees (default: %(default)s)')

    group = parser.add_argument_group('Training parameters')
    group.add_argument('-n', '--num-epochs', type=int, default=10, help='Number of training epochs (default: %(default)s)')
    group.add_argument('-b','--batch-size', type=int, default=100, help='Minibatch size (default: %(default)s)')
    group.add_argument('--wd', type=float, default=0, help='Weight decay in Adam optimizer (default: %(default)s)')
    group.add_argument('--lr', type=float, default=1e-3, help='Learning rate in Adam optimizer (default: %(default)s)')
    group.add_argument('--beta', default=1.0, help='Choice of beta schedule or a constant for KLD weight (default: %(default)s)')
    group.add_argument('--beta-control', type=float, help='KL-Controlled VAE gamma. Beta is KL target. (default: %(default)s)')
    group.add_argument('--equivariance', type=float, help='Strength of equivariance loss (default: %(default)s)')
    group.add_argument('--equivariance-end-it', type=int, default=100000, help='It at which equivariance max (default: %(default)s)')
    group.add_argument('--l-start', type=int,default=12, help='Starting L radius (default: %(default)s)')
    group.add_argument('--l-end', type=int, default=20, help='End L radius (default: %(default)s)')
    group.add_argument('--l-end-it',type=int,default=100000, help='default: %(default)s')

    group = parser.add_argument_group('Encoder Network')
    group.add_argument('--qlayers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--qdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--encode-mode', default='resid', choices=('conv','resid','mlp','tilt'), help='Type of encoder network (default: %(default)s)')
    group.add_argument('--zdim', type=int, default=1, help='Dimension of latent variable')

    group = parser.add_argument_group('Decoder Network')
    group.add_argument('--players', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--pdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    return parser

def eval_volume(model, lattice, D, zval, rnorm):
    '''Evaluate the model on a nz x ny x nx lattice'''
    zdim = len(zval)
    z = torch.zeros(D**2,zdim, dtype=torch.float32)
    z += torch.tensor(zval, dtype=torch.float32)

    vol_f = np.zeros((D,D,D),dtype=np.float32)
    assert not model.training
    # evaluate the volume by zslice to avoid memory overflows
    for i, dz in enumerate(np.linspace(-1,1,D,endpoint=False)):
        x = lattice.coords + torch.tensor([0,0,dz], dtype=torch.float32)
        x = torch.cat((x,z),dim=-1)
        with torch.no_grad():
            y = model.decoder.decode(x)
            y = y[...,0] - y[...,1]
            #y = model.decoder(x)
            y = y.view(D,D).cpu().numpy()
        vol_f[i] = y
    vol = fft.ihtn_center(vol_f*rnorm[1]+rnorm[0])
    return vol, vol_f

def main(args):
    log(args)
    t1 = dt.now()
    if args.outdir is not None and not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # set the random seed
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    ## set the device
    use_cuda = torch.cuda.is_available()
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
    particles, _, _ = mrc.parse_mrc(args.particles)
    Nimg, ny, nx = particles.shape
    assert ny == nx
    nz = max(nx,ny)
    log('Loaded {} {}x{} images'.format(Nimg, ny, nx))
    particles = np.asarray([fft.ht2_center(img).astype(np.float32) for img in particles])
    assert particles.shape == (Nimg,ny,nx)
    rnorm  = [np.mean(particles), np.std(particles)]
    log('Particle FT stack mean, std: {} +/- {}'.format(*rnorm))
    rnorm[0] = 0
    log('Normalizing FT by mean, std: {} +/- {}'.format(*rnorm))
    particles = (particles - rnorm[0])/rnorm[1]

    # load particles from tilt series
    if args.tilt is not None:
        particles_tilt, _, _ = mrc.parse_mrc(args.tilt)
        assert particles_tilt.shape == (Nimg, ny, nx), 'Tilt series pair must have same dimensions as untilted particles'
        log('Loaded {} {}x{} tilt series images'.format(Nimg, ny, nx))
        particles_tilt = np.asarray([fft.ht2_center(img).astype(np.float32) for img in particles_tilt])
        particles_tilt = (particles_tilt - rnorm[0])/rnorm[1]
    
        theta = args.tilt_deg*np.pi/180
        tilt = np.array([[1.,0.,0.],
                        [0, np.cos(theta), -np.sin(theta)],
                        [0, np.sin(theta), np.cos(theta)]]).astype(np.float32)
        tilt = torch.tensor(tilt)
        assert args.encode_mode == 'tilt'
    else:
        tilt = None

    lattice = Lattice(nx)
    model = HetVAE(lattice, nx*ny, args.qlayers, args.qdim, args.players, args.pdim,
                args.zdim, encode_mode=args.encode_mode)
    bnb = BNBHetOpt(model, lattice, tilt)

    if args.equivariance:
        assert args.equivariance > 0, 'Regularization weight must be positive'
        equivariance_lambda = LinearSchedule(0, args.equivariance, 10000, args.equivariance_end_it)
        equivariance_loss = EquivarianceLoss(model,ny,nx)

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
    for epoch in range(start_epoch, num_epochs):
        kld_accum = 0
        gen_loss_accum = 0
        loss_accum = 0
        eq_loss_accum = 0
        batch_it = 0 
        num_batches = np.ceil(Nimg / args.batch_size).astype(int)
        z_accum = np.zeros(args.zdim)
        bnb_result = []
        for minibatch_i in np.array_split(np.random.permutation(Nimg),num_batches):
            batch_it += len(minibatch_i)
            global_it = Nimg*epoch + batch_it

            y = torch.from_numpy(particles[minibatch_i])
            if use_cuda: y = y.cuda()
            if tilt is not None:
                yt = torch.from_numpy(particles_tilt[minibatch_i])
                if use_cuda: yt = yt.cuda()
            else:
                yt = None

            # predict encoding
            #input_ = torch.stack((y, yt),1) if tilt is not None else y
            input_ = (y,yt) if tilt is not None else y
            mu, logvar = model.encode(input_)
            z = model.reparameterize(mu, logvar)

            # equivariance loss
            if args.equivariance:
                lamb = equivariance_lambda(global_it)
                eq_loss = equivariance_loss(y, mu)
            else:
                lamb, eq_loss = 0, 0 

            # find the optimal orientation for each image
            L = Lsched(global_it)
            if L: L = int(L)
            model.eval()
            with torch.no_grad():
                rot = bnb.opt_theta(y, z, yt, L=L)
            model.train()

            # train the decoder
            y_recon = model(rot, z)
            y_recon = y_recon.view(-1, ny, nx)
            gen_loss = F.mse_loss(y_recon,y)
            if tilt is not None:
                y_recon_tilt = model(bnb.tilt @ rot, z)
                y_recon_tilt = y_recon_tilt.view(-1, ny, nx)
                gen_loss = .5*gen_loss + F.mse_loss(y_recon_tilt,yt)

            kld = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())

            beta = beta_schedule(global_it)
            if args.beta_control is None:
                loss = gen_loss + beta*kld/(nx*ny)
            else:
                loss = gen_loss + args.beta_control*(beta-kld)**2/(nx*ny)

            if args.equivariance:
                loss += lamb*eq_loss

            if torch.isnan(kld):
                log(mu[0])
                log(logvar[0])
                raise RuntimeError('KLD is nan')


            loss.backward()
            optim.step()
            optim.zero_grad()
            
            # logging
            bnb_result.append((minibatch_i,rot.detach().cpu().numpy()))
            z_accum += mu.detach().cpu().numpy().mean()*len(minibatch_i)
            kld_accum += kld.item()*len(minibatch_i)
            gen_loss_accum += gen_loss.item()*len(minibatch_i)
            loss_accum += loss.item()*len(minibatch_i)
            if args.equivariance:eq_loss_accum += eq_loss.item()*len(minibatch_i)

            if batch_it % args.log_interval == 0:
                eq_log = 'equivariance={:.4f}, lambda={:.4f}, '.format(eq_loss.item(), lamb) if args.equivariance else ''
                log('# [Train Epoch: {}/{}] [{}/{} images] gen loss={:.4f}, kld={:.4f}, beta={:.4f}, {}loss={:.4f}'.format(epoch+1, num_epochs, batch_it, Nimg, gen_loss.item(), kld.item(), beta, eq_log, loss.item()))
        eq_log = 'equivariance = {:.4f}, '.format(eq_loss_accum/Nimg) if args.equivariance else ''
        log('# =====> Epoch: {} Average gen loss = {:.4}, KLD = {:.4f}, {}total loss = {:.4f}'.format(epoch+1, gen_loss_accum/Nimg, kld_accum/Nimg, eq_log, loss_accum/Nimg))


        if args.checkpoint and epoch % args.checkpoint == 0:
            model.eval()
            vol, vol_f = eval_volume(model, lattice, nx, z_accum/Nimg, rnorm)
            mrc.write('{}/reconstruct.{}.mrc'.format(args.outdir,epoch), vol.astype(np.float32))
            path = '{}/weights.{}.pkl'.format(args.outdir,epoch)
            torch.save({
                'norm':rnorm,
                'epoch':epoch,
                'model_state_dict':model.state_dict(),
                'optimizer_state_dict':optim.state_dict(),
                'bnb_result':bnb_result
                }, path)
            model.train()

    ## save model weights and evaluate the model on 3D lattice
    model.eval()
    vol, vol_f = eval_volume(model, lattice, nx, z_accum/Nimg, rnorm)
    mrc.write('{}/reconstruct.mrc'.format(args.outdir), vol.astype(np.float32))
    path = '{}/weights.pkl'.format(args.outdir)
    torch.save({
        'norm':rnorm,
        'epoch':epoch,
        'model_state_dict':model.state_dict(),
        'optimizer_state_dict':optim.state_dict(),
        'bnb_result':bnb_result
        }, path)
    
    td = dt.now()-t1
    log('Finsihed in {} ({} per epoch)'.format(td, td/num_epochs))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

