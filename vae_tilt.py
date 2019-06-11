'''
VAE-based homogeneous reconstruction with a tilt series
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
import lie_tools
import dataset

from lattice import Lattice
from models import TiltVAE
from beta_schedule import get_beta_schedule, LinearSchedule

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrcs)')
    parser.add_argument('particles_tilt', help='Particle stack file for tilt pair (.mrcs)')
    parser.add_argument('--tilt', type=float, default=45, help='Right-handed x-axis tilt offset in degrees (default: %(default)s)')
    parser.add_argument('-o', '--outdir', type=os.path.abspath, required=True, help='Output directory to save model')
    parser.add_argument('--norm', type=float, nargs=2, default=None, help='Data normalization as shift, 1/scale (default: mean, std of dataset)')
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
    group.add_argument('--beta', default=1.0, help='Choice of beta schedule or a constant for KLD weight (default: %(default)s)')
    group.add_argument('--beta-control', type=float, help='KL-Controlled VAE gamma. Beta is KL target. (default: %(default)s)')
    group.add_argument('--no-trans', action='store_true', help="Don't translate images")
       
    group = parser.add_argument_group('Encoder Network')
    group.add_argument('--qlayers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--qdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--enc-mask', type=int, help='Radius of circular mask for image encoder')

    group = parser.add_argument_group('Decoder Network')
    group.add_argument('--players', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--pdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    return parser

def loss_function(recon_y, recon_y_tilt, y, yt, w_eps, z_std, t_mu, t_logvar):
    gen_loss = F.mse_loss(recon_y, y)*.5 + F.mse_loss(recon_y_tilt, yt)*.5   
    cross_entropy = torch.tensor([np.log(8*np.pi**2)], device=y.device) # cross entropy between gaussian and uniform on SO3
    entropy = lie_tools.so3_entropy(w_eps,z_std)
    kld1 = (cross_entropy - entropy).mean()
    if t_mu is not None:
        kld2 = -0.5 * torch.mean(1 + t_logvar - t_mu.pow(2) - t_logvar.exp())
    else: kld2 = 0.0
    #assert kld > 0
    return gen_loss, kld1 + kld2

def train(model, optim, D, minibatch, beta, beta_control=None):
    y, yt = minibatch
    model.train()
    optim.zero_grad()
    # train the model
    y_recon, y_recon_tilt, y, yt, z_mu, z_std, w_eps, t_mu, t_logvar = model(y, yt) 
    gen_loss, kld = loss_function(y_recon, y_recon_tilt, y, yt, w_eps, z_std, t_mu, t_logvar)
    if torch.isnan(kld):
        log(w_eps[0])
        log(z_std[0])
        raise RuntimeError('KLD is nan')
    if args.beta_control is None:
        loss = gen_loss + beta*kld/D**2
    else:
        loss = gen_loss + args.beta_control*(beta-kld)**2/D**2
    loss.backward()
    optim.step()
    return gen_loss.item(), kld.item(), loss.item()

def save_checkpoint(model, optim, D, epoch, norm, out_mrc, out_weights):
    model.eval()
    vol = model.eval_volume(norm)
    mrc.write(out_mrc, vol.astype(np.float32))
    torch.save({
        'epoch':epoch,
        'model_state_dict':model.state_dict(),
        'optimizer_state_dict':optim.state_dict(),
        }, out_weights)

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
    data = dataset.TiltMRCData(args.particles, args.particles_tilt, norm=args.norm)
    Nimg = data.N
    D = data.D

    theta = args.tilt*np.pi/180
    tilt = np.array([[1.,0.,0.],
                    [0, np.cos(theta), -np.sin(theta)],
                    [0, np.sin(theta), np.cos(theta)]]).astype(np.float32)

    lattice = Lattice(D)
    if args.enc_mask: args.enc_mask = lattice.get_circular_mask(args.enc_mask)
    model = TiltVAE(lattice, tilt, args.qlayers, args.qdim, args.players, args.pdim, no_trans=args.no_trans, enc_mask=args.enc_mask)
    log(model)
    log('{} parameters in model'.format(sum(p.numel() for p in model.parameters() if p.requires_grad)))


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
    data_generator = DataLoader(data, batch_size=args.batch_size, shuffle=True)
    num_epochs = args.num_epochs
    for epoch in range(start_epoch, num_epochs):
        gen_loss_accum = 0
        loss_accum = 0
        kld_accum = 0
        batch_it = 0 
        for minibatch in data_generator:
            minibatch = (minibatch[0].to(device), minibatch[1].to(device))
            batch_it += len(minibatch[0])
            global_it = Nimg*epoch + batch_it

            beta = beta_schedule(global_it)

            gen_loss, kld, loss = train(model, optim, D, minibatch, beta, args.beta_control)
            # logging
            gen_loss_accum += gen_loss*len(minibatch[0])
            kld_accum += kld*len(minibatch[0])
            loss_accum += loss*len(minibatch[0])

            if batch_it % args.log_interval == 0:
                log('# [Train Epoch: {}/{}] [{}/{} images] gen loss={:.4f}, kld={:.4f}, beta={:.4f}, loss={:.4f}'.format(epoch+1, num_epochs, batch_it, Nimg, gen_loss, kld, beta, loss))
        log('# =====> Epoch: {} Average gen loss = {:.4}, KLD = {:.4f}, total loss = {:.4f}'.format(epoch+1, gen_loss_accum/Nimg, kld_accum/Nimg, loss_accum/Nimg))

        if args.checkpoint and epoch % args.checkpoint == 0:
            out_mrc = '{}/reconstruct.{}.mrc'.format(args.outdir,epoch)
            out_weights = '{}/weights.{}.pkl'.format(args.outdir,epoch)
            save_checkpoint(model, optim, D, epoch, data.norm, out_mrc, out_weights)

    ## save model weights and evaluate the model on 3D lattice
    out_mrc = '{}/reconstruct.mrc'.format(args.outdir)
    out_weights = '{}/weights.pkl'.format(args.outdir)
    save_checkpoint(model, optim, D, epoch, data.norm, out_mrc, out_weights)
    
    td = dt.now()-t1
    log('Finsihed in {} ({} per epoch)'.format(td, td/num_epochs))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

