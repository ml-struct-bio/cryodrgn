'''
Homogeneous NN reconstruction with hierarchical pose optimization
'''

import numpy as np
import sys, os
import argparse
import pickle
import glob
from datetime import datetime as dt

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader

from cryodrgn import mrc
from cryodrgn import utils
from cryodrgn import fft
from cryodrgn import dataset
from cryodrgn import lie_tools
from cryodrgn import ctf

from cryodrgn.lattice import Lattice
from cryodrgn.pose_search import PoseSearch
from cryodrgn import models

log = utils.log
vlog = utils.vlog

def debug_signal_handler(signal, frame):
    import pdb
    pdb.set_trace()
import signal
signal.signal(signal.SIGINT, debug_signal_handler)


def add_args(parser):
    parser.add_argument('particles', type=os.path.abspath, help='Particle stack file (.mrcs)')
    parser.add_argument('-o', '--outdir', type=os.path.abspath, required=True, help='Output directory to save model')
    parser.add_argument('--ctf', metavar='pkl', type=os.path.abspath, help='CTF parameters (.pkl)')
    parser.add_argument('--norm', type=float, nargs=2, default=None, help='Data normalization as shift, 1/scale (default: mean, std of dataset)')
    parser.add_argument('--load', help='Initialize training from a checkpoint')
    parser.add_argument('--load-poses', type=os.path.abspath, help='Initialize training from a checkpoint')
    parser.add_argument('--checkpoint', type=int, default=1, help='Checkpointing interval in N_EPOCHS (default: %(default)s)')
    parser.add_argument('--log-interval', type=int, default=1000, help='Logging interval in N_IMGS (default: %(default)s)')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')
    parser.add_argument('--seed', type=int, default=np.random.randint(0,100000), help='Random seed')
    parser.add_argument('--uninvert-data', dest='invert_data', action='store_false', help='Do not invert data sign')
    parser.add_argument('--window', action='store_true', help='Real space windowing of dataset')
    parser.add_argument('--window-r', type=float, default=.85,  help='Windowing radius (default: %(default)s)')
    parser.add_argument('--ind', type=os.path.abspath, help='Filter indices')

    group = parser.add_argument_group('Tilt series')
    group.add_argument('--tilt', help='Particle stack file (.mrcs)')
    group.add_argument('--tilt-deg', type=float, default=45, help='X-axis tilt offset in degrees (default: %(default)s)')

    group = parser.add_argument_group('Training parameters')
    group.add_argument('--t-extent', type=float, default=10, help='+/- pixels to search over translations')
    group.add_argument('--t-ngrid', type=float, default=7, help='Initial grid size for translations')
    group.add_argument('--t-xshift', type=float, default=0)
    group.add_argument('--t-yshift', type=float, default=0)
    group.add_argument('--no-trans', action='store_true', help="Don't search over translations")
    group.add_argument('--pretrain', type=int, default=10000, help='Number of initial iterations with random poses (default: %(default)s)')
    group.add_argument('--ps-freq', type=int, default=5, help='Frequency of pose inference (default: every %(default)s epochs)')
    group.add_argument('-n', '--num-epochs', type=int, default=30, help='Number of training epochs (default: %(default)s)')
    group.add_argument('-b','--batch-size', type=int, default=10, help='Minibatch size (default: %(default)s)')
    group.add_argument('--wd', type=float, default=0, help='Weight decay in Adam optimizer (default: %(default)s)')
    group.add_argument('--lr', type=float, default=1e-4, help='Learning rate in Adam optimizer (default: %(default)s)')
    group.add_argument('--reset-model-every', type=int, help="If set, reset the model every N epochs")
    group.add_argument('--reset-optim-every', type=int, help="If set, reset the optimizer every N epochs")
    group.add_argument('--reset-optim-after-pretrain', type=int, help="If set, reset the optimizer every N epochs")

    group = parser.add_argument_group('Pose search parameters')
    group.add_argument('--l-start', type=int,default=12, help='Starting L radius (default: %(default)s)')
    group.add_argument('--l-end', type=int, default=32, help='End L radius (default: %(default)s)')
    group.add_argument('--niter', type=int, default=4, help='Number of iterations of grid subdivision')
    group.add_argument('--l-ramp-epochs',type=int,default=25, help='Number of epochs to ramp up to --l-end (default: %(default)s)')
    group.add_argument('--probabilistic', action='store_true', help='Use probabilistic bound')
    group.add_argument('--nkeptposes', type=int, default=8, help="Number of poses to keep at each refinement interation during branch and bound")
    group.add_argument('--base-healpy', type=int, default=2, help="Base healpy grid for pose search. Higher means exponentially higher resolution.")
    group.add_argument("--pose-model-update-freq", type=int, help="If set, only update the model used for pose search every N examples.")

    group = parser.add_argument_group('Network Architecture')
    group.add_argument('--layers', type=int, default=3, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--dim', type=int, default=256, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--l-extent', type=float, default=0.5, help='Coordinate lattice size (if not using positional encoding) (default: %(default)s)')
    group.add_argument('--pe-type', choices=('geom_ft','geom_full','geom_lowf','geom_nohighf','linear_lowf','gaussian','none'), default='gaussian', help='Type of positional encoding (default: %(default)s)')
    group.add_argument('--pe-dim', type=int, help='Num sinusoid features in positional encoding (default: D/2)')
    group.add_argument('--domain', choices=('hartley','fourier'), default='fourier', help='Volume decoder representation (default: %(default)s)')
    group.add_argument('--activation', choices=('relu','leaky_relu'), default='relu', help='Activation (default: %(default)s)')
    group.add_argument('--feat-sigma', type=float, default=0.5, help="Scale for random Gaussian features")

    return parser

def save_checkpoint(model, lattice, pose, optim, epoch, norm, out_mrc, out_weights, out_poses):
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
        pickle.dump(pose, f)

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

def sort_poses(pose):
    ind = [x[0] for x in pose]
    ind = np.concatenate(ind)
    rot = [x[1][0] for x in pose]
    rot = np.concatenate(rot)

    rot = rot[np.argsort(ind)]
    if len(pose[0][1]) == 2:
        trans = [x[1][1] for x in pose]
        trans = np.concatenate(trans)
        trans = trans[np.argsort(ind)]
        return (rot,trans)
    return (rot,)


def sort_base_poses(pose):
    ind, data = zip(*pose)
    ind = np.concatenate(ind)
    data = np.concatenate(data)
    return data[np.argsort(ind)]


def train(model, lattice, ps, optim, batch, tilt_rot=None, no_trans=False, poses=None, base_pose=None, ctf_params=None):
    y, yt = batch
    B = y.size(0)
    D = lattice.D

    if ctf_params is not None:
        freqs = lattice.freqs2d.unsqueeze(0).expand(B,*lattice.freqs2d.shape)/ctf_params[:,0].view(B,1,1)
        ctf_i = ctf.compute_ctf(freqs, *torch.split(ctf_params[:,1:], 1, 1)).view(B,D,D)
    else: ctf_i = None

    # pose inference
    if poses is not None:
        rot = poses[0]
        if not no_trans: trans = poses[1]
    else: # BNB
        model.eval()
        with torch.no_grad():
            rot, trans, base_pose = ps.opt_theta_trans(y, images_tilt=yt, init_poses=base_pose, ctf_i=ctf_i)
            base_pose = base_pose.detach().cpu().numpy()

    # reconstruct circle of pixels instead of whole image
    mask = lattice.get_circular_mask(lattice.D//2)
    # mask = lattice.get_circular_mask(ps.Lmax)
    def gen_slice(R):
        slice_ = model(lattice.coords[mask] @ R).view(B,-1)
        if ctf_params is not None:
            slice_ *= ctf_i.view(B,-1)[:,mask]
        return slice_
    def translate(img):
        img = lattice.translate_ht(img, trans.unsqueeze(1), mask)
        return img.view(B,-1)

    # Train model
    model.train()
    optim.zero_grad()

    y = y.view(B,-1)[:, mask]
    if tilt_rot is not None: yt = yt.view(B,-1)[:, mask]
    if not no_trans:
        y = translate(y)
        if tilt_rot is not None: yt = translate(yt)

    if tilt_rot is not None:
        loss = .5*F.mse_loss(gen_slice(rot), y) + .5*F.mse_loss(gen_slice(tilt_rot @ rot), yt)
    else:
        loss = F.mse_loss(gen_slice(rot), y)
    loss.backward()
    optim.step()
    save_pose = [rot.detach().cpu().numpy()]
    if not no_trans:
        save_pose.append(trans.detach().cpu().numpy())
    return loss.item(), save_pose, base_pose

def get_latest(args, flog):
    # assumes args.num_epochs > latest checkpoint
    flog('Detecting latest checkpoint...') 
    weights = [f'{args.outdir}/weights.{i}.pkl' for i in range(args.num_epochs)]
    weights = [f for f in weights if os.path.exists(f)]
    args.load = weights[-1]
    flog(f'Loading {args.load}')
    i = args.load.split('.')[-2]
    args.load_poses = f'{args.outdir}/pose.{i}.pkl'
    assert os.path.exists(args.load_poses) # Might need to relax this assert
    flog(f'Loading {args.load_poses}')
    return args

def make_model(args, D):
    activation={"relu": nn.ReLU, "leaky_relu": nn.LeakyReLU}[args.activation]
    return models.get_decoder(3, D, args.layers, args.dim, args.domain, args.pe_type, enc_dim=args.pe_dim, activation=activation, feat_sigma=args.feat_sigma)

def main(args):
    t1 = dt.now()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    LOG = f'{args.outdir}/run.log'
    def flog(msg): # HACK: switch to logging module
        return utils.flog(msg, LOG)
    if args.load == 'latest':
        args = get_latest(args, flog)
    flog(' '.join(sys.argv))
    flog(args)

    # set the random seed
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    ## set the device
    use_cuda = torch.cuda.is_available()
    device = torch.device('cuda' if use_cuda else 'cpu')
    flog('Use cuda {}'.format(use_cuda))
    if not use_cuda:
        flog('WARNING: No GPUs detected')

    # load the particles
    if args.ind is not None: 
        flog('Filtering image dataset with {}'.format(args.ind))
        args.ind = pickle.load(open(args.ind,'rb'))
    if args.tilt is None:
        data = dataset.MRCData(args.particles, norm=args.norm, invert_data=args.invert_data, ind=args.ind, window=args.window, window_r=args.window_r, flog=flog)
        tilt = None
    else:
        data = dataset.TiltMRCData(args.particles, args.tilt, norm=args.norm, invert_data=args.invert_data, ind=args.ind, window=args.window, window_r=args.window_r, flog=flog)
        tilt = torch.tensor(utils.xrot(args.tilt_deg).astype(np.float32), device=device)
    D = data.D
    Nimg = data.N

    # load ctf
    if args.ctf is not None:
        flog('Loading ctf params from {}'.format(args.ctf))
        ctf_params = ctf.load_ctf_for_training(D-1, args.ctf)
        if args.ind is not None: ctf_params = ctf_params[args.ind] 
        assert ctf_params.shape == (Nimg, 8)
        ctf_params = torch.tensor(ctf_params, device=device)
    else: ctf_params = None

    # instantiate model
    lattice = Lattice(D, extent=0.5, device=device)
    model = make_model(args, D)
    model.to(device)

    if args.pose_model_update_freq:
        pose_model = make_model(args, D)
        pose_model.to(device)
        pose_model.eval()
    else:
        pose_model = model

    if args.no_trans:
        raise NotImplementedError()
    else:
        ps = PoseSearch(pose_model, lattice, args.l_start, args.l_end, tilt,
                        t_extent=args.t_extent, t_ngrid=args.t_ngrid, niter=args.niter,
                        nkeptposes=args.nkeptposes, base_healpy=args.base_healpy,
                        t_xshift=args.t_xshift, t_yshift=args.t_yshift, device=device)
    flog(model)
    flog('{} parameters in model'.format(sum(p.numel() for p in model.parameters() if p.requires_grad)))
    optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

    if args.load:
        args.pretrain = 0
        flog('Loading checkpoint from {}'.format(args.load))
        checkpoint = torch.load(args.load)
        model.load_state_dict(checkpoint['model_state_dict'])
        optim.load_state_dict(checkpoint['optimizer_state_dict'])
        start_epoch = checkpoint['epoch']+1
        assert args.num_epochs > start_epoch
        model.train()
        if args.load_poses:
            sorted_poses = utils.load_pkl(args.load_poses)
            sorted_base_poses = None # TODO: need to save base_poses if we are going to use it
    else:
        start_epoch = 0

    data_iterator = DataLoader(data, batch_size=args.batch_size, shuffle=True)

    # pretrain decoder with random poses
    global_it = 0
    pretrain_epoch = 0
    flog('Using random poses for {} iterations'.format(args.pretrain))
    while global_it < args.pretrain:
        for batch in data_iterator:
            global_it += len(batch[0])
            batch = (batch[0].to(device), None) if tilt is None else (batch[0].to(device), batch[1].to(device))
            loss = pretrain(model, lattice, optim, batch, tilt=ps.tilt)
            if global_it % args.log_interval == 0:
                flog(f'[Pretrain Iteration {global_it}] loss={loss:4f}')
            if global_it > args.pretrain:
                break
    out_mrc = '{}/pretrain.reconstruct.mrc'.format(args.outdir)
    model.eval()
    vol = model.eval_volume(lattice.coords, lattice.D, lattice.extent, data.norm)
    mrc.write(out_mrc, vol.astype(np.float32))

    # reset model after pretraining
    if args.reset_optim_after_pretrain:
        flog(">> Resetting optim after pretrain")
        optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

    # training loop
    cc = 0
    if args.pose_model_update_freq:
        pose_model.load_state_dict(model.state_dict())
    for epoch in range(start_epoch, args.num_epochs):
        t2 = dt.now()
        batch_it = 0
        loss_accum = 0
        poses, base_poses = [], []

        if args.l_ramp_epochs > 0:
            Lramp = args.l_start + int(epoch / args.l_ramp_epochs * (args.l_end - args.l_start))
            ps.Lmin = min(Lramp, args.l_start)
            ps.Lmax = min(Lramp, args.l_end)

        if args.reset_model_every and (epoch - 1) % args.reset_model_every == 0:
            flog(">> Resetting model")
            model = make_model(args, D)

        if args.reset_optim_every and (epoch - 1) % args.reset_optim_every == 0:
            flog(">> Resetting optim")
            optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

        if epoch % args.ps_freq != 0:
            flog('Using previous iteration poses')
        for batch in data_iterator:
            ind = batch[-1]
            ind_np = ind.cpu().numpy()
            batch = (batch[0].to(device), None) if tilt is None else (batch[0].to(device), batch[1].to(device))
            batch_it += len(batch[0])

            # train the model
            if epoch % args.ps_freq != 0:
                p = [torch.tensor(x[ind_np], device=device) for x in sorted_poses]
                #bp = sorted_base_poses[ind_np]
                bp = None
            else:
                p, bp = None, None

            cc += len(batch[0])
            if args.pose_model_update_freq and cc > args.pose_model_update_freq:
                pose_model.load_state_dict(model.state_dict())
                cc = 0
                
            c = ctf_params[ind] if ctf_params is not None else None
            loss_item, pose, base_pose = train(model, lattice, ps, optim, batch, tilt, args.no_trans, poses=p, base_pose=bp, ctf_params=c)
            poses.append((ind.cpu().numpy(),pose))
            base_poses.append((ind_np, base_pose))
            # logging
            loss_accum += loss_item*len(batch[0])
            if batch_it % args.log_interval == 0:
                log('# [Train Epoch: {}/{}] [{}/{} images] loss={:.4f}'.format(epoch+1, args.num_epochs, batch_it, Nimg, loss_item))
        flog('# =====> Epoch: {} Average loss = {:.4}; Finished in {}'.format(epoch+1, loss_accum/Nimg, dt.now() - t2))

        # sort pose
        sorted_poses = sort_poses(poses) if poses else None
        #sorted_base_poses = sort_base_poses(base_poses)

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
    flog('Finished in {} ({} per epoch)'.format(td, td/(args.num_epochs-start_epoch)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    utils._verbose = args.verbose
    main(args)

