'''
'''
import numpy as np
import sys, os
import argparse
import pickle
from datetime import datetime as dt
import matplotlib.pyplot as plt 
import pprint
import shutil

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset
from sklearn.model_selection import train_test_split

import cryodrgn
from cryodrgn import mrc
from cryodrgn import utils
from cryodrgn import fft
from cryodrgn import lie_tools
from cryodrgn import config
from cryodrgn.lattice import Lattice
from cryodrgn.models import HetOnlyVAE, ResidLinearMLP

log = utils.log
vlog = utils.vlog

def add_args(parser):
    parser.add_argument('workdir', type=os.path.abspath, help='Directory with cryoDRGN results')
    parser.add_argument('epoch', type=int, help='Epoch number N to analyze (0-based indexing, corresponding to z.N.pkl, weights.N.pkl)')
    parser.add_argument('--device', type=int, help='Optionally specify CUDA device')
    parser.add_argument('--landscape-dir', type=os.path.abspath, help='Landscape analysis directory (default: [workdir]/landscape.[epoch])')
    parser.add_argument('-o','--outdir', type=os.path.abspath, help='Output directory (default: [workdir]/landscape.[epoch]/landscape_full)')
    parser.add_argument('--seed', default=0, type=int, help='Random seed (default: %(default)s)')

    group = parser.add_argument_group('Volume generation arguments')
    group.add_argument('-N', type=int, default=10000, help='Number of training volumes to generate (default: %(default)s)')
    group.add_argument('--flip', action='store_true', help='Flip handedness')
    group.add_argument('-d','--downsample', type=int, default=128, help='Downsample volumes to this box size (pixels) (default: %(default)s)')
    group.add_argument('--skip-vol', action='store_true', help='Skip generation of volumes')


    group = parser.add_argument_group('Volume mapping arguments')
    group.add_argument('--batch-size', type=int, default=64, metavar='N',
                        help='input batch size for training (default: 64)')
    group.add_argument('--test-batch-size', type=int, default=1000, metavar='N',
                        help='input batch size for testing (default: 1000)')
    group.add_argument('--epochs', type=int, default=50, metavar='N',
                        help='number of epochs to train (default: 14)')
    group.add_argument('--lr', type=float, default=1e-3, metavar='LR',
                        help='learning rate (default: 1e-3)')
    group.add_argument('--dim', type=int, default=512, metavar='N',
                        help='MLP hidden layer dimension (default: 512)')
    group.add_argument('--layers', type=int, default=3, metavar='N',
                        help='MLP number of hidden layers (default: 3)')
    return parser

def train(args, model, device, train_loader, optimizer, epoch):
    model.train()
    for batch_idx, (data, target) in enumerate(train_loader):
        data, target = data.to(device), target.to(device)
        optimizer.zero_grad()
        output = model(data)
        loss = F.mse_loss(output, target)
        loss.backward()
        optimizer.step()
        if batch_idx % 10 == 0:
            print('Train Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(
                epoch, batch_idx * len(data), len(train_loader.dataset),
                100. * batch_idx / len(train_loader), loss.item()))


def test(model, device, test_loader):
    model.eval()
    test_loss = 0
    with torch.no_grad():
        for data, target in test_loader:
            data, target = data.to(device), target.to(device)
            output = model(data)
            test_loss += len(data)*F.mse_loss(output, target).item()  # sum up batch loss

    test_loss /= len(test_loader.dataset)

    print('\nTest set: Average loss: {:.4f}\n'.format(test_loss))

class MyDataset(Dataset):
    def __init__(self, x, y):
        assert len(x) == len(y)
        self.x = x
        self.y = y

    def __len__(self):
        return len(self.x)

    def __getitem__(self, idx):
        return self.x[idx], self.y[idx]

def generate_and_map_volumes(zfile, cfg_pkl, weights, mask_mrc, pca_obj_pkl, landscape_dir, outdir, args):
    # Sample z
    log(f'Sampling {args.N} particles from {zfile}')
    np.random.seed(args.seed)
    z_all = utils.load_pkl(zfile)
    ind = np.array(sorted(np.random.choice(len(z_all), args.N, replace=False)))
    z_sample = z_all[ind]
    utils.save_pkl(z_sample, f'{outdir}/z.sampled.pkl')
    utils.save_pkl(ind, f'{outdir}/ind.sampled.pkl')
    log(f'Saved {outdir}/z.sampled.pkl')

    # Set the device
    assert torch.cuda.is_available(), "No GPUs detected"
    torch.set_default_tensor_type(torch.cuda.FloatTensor)

    cfg = config.update_config_v1(cfg_pkl)
    log('Loaded configuration:')
    pprint.pprint(cfg)

    D = cfg['lattice_args']['D'] # image size + 1
    norm = cfg['dataset_args']['norm']

    # Load landscape analysis inputs
    mask, _ = mrc.parse_mrc(mask_mrc)
    mask = mask.astype(bool)
    if args.downsample: 
        assert mask.shape == (args.downsample,)*3
    else:
        assert mask.shape == (D-1,D-1,D-1)
    log(f'{mask.sum()} voxels in the mask')

    pca = utils.load_pkl(pca_obj_pkl)

    # Load model weights
    log('Loading weights from {}'.format(weights))
    model, lattice = HetOnlyVAE.load(cfg, weights)
    model.eval()

    # Set z 
    z = z_sample.astype(np.float32)
   
    # Generate volumes
    log(f'Generating {len(z)} volume embeddings')
    t1 = dt.now()
    embeddings = []
    for i,zz in enumerate(z):
        if i % 100 == 0:
            log(i)
        if args.downsample:
            extent = lattice.extent * (args.downsample/(D-1))
            vol = model.decoder.eval_volume(lattice.get_downsample_coords(args.downsample+1), 
                                            args.downsample+1, extent, norm, zz) 
        else:
            vol = model.decoder.eval_volume(lattice.coords, lattice.D, lattice.extent, norm, zz) 
        if args.flip:
            vol = vol[::-1]
        embeddings.append(pca.transform(vol[mask].reshape(1,-1)))

    embeddings = np.array(embeddings).reshape(len(z),-1).astype(np.float32)

    td = dt.now()-t1
    log(f'Finished generating {args.N} volumes in {td}')
    return z, embeddings

def train_model(x, y, outdir, zfile, args):
    use_cuda = torch.cuda.is_available()
    torch.manual_seed(args.seed)
    device = torch.device("cuda" if use_cuda else "cpu")

    train_kwargs = {'batch_size': args.batch_size}
    test_kwargs = {'batch_size': args.test_batch_size}
    if use_cuda:
        cuda_kwargs = {'num_workers': 1,
                       'pin_memory': True,
                       'shuffle': True}
        train_kwargs.update(cuda_kwargs)
        test_kwargs.update(cuda_kwargs)

    # Load dataset
    x_train, x_test, y_train, y_test = train_test_split(x,y,test_size=.25, random_state=args.seed)

    train_dataset = MyDataset(x_train, y_train)
    test_dataset = MyDataset(x_test, y_test)

    train_loader = torch.utils.data.DataLoader(train_dataset,**train_kwargs)
    test_loader = torch.utils.data.DataLoader(test_dataset, **test_kwargs)

    model = ResidLinearMLP(x.shape[1], args.layers, args.dim, y.shape[1], nn.ReLU).to(device)
    log(model)
    optimizer = optim.Adam(model.parameters(), lr=args.lr)

    # Train
    for epoch in range(1, args.epochs+1):
        train(args, model, device, train_loader, optimizer, epoch)
        test(model, device, test_loader)

    # Evaluate
    model.eval()
    yhat_all = []
    eval_dataset = utils.load_pkl(zfile).astype(np.float32)
    with torch.no_grad():
        for x in np.array_split(eval_dataset, args.test_batch_size):
            x = torch.tensor(x,device=device)
            yhat = model(x)
            yhat_all.append(yhat.detach().cpu().numpy())
    yhat_all = np.concatenate(yhat_all)
    torch.save(model.state_dict(), f"{outdir}/model.pt")
    return yhat_all


def main(args):
    t1 = dt.now()
    log(args)
    E = args.epoch
    workdir = args.workdir
    zfile = f'{workdir}/z.{E}.pkl'
    weights = f'{workdir}/weights.{E}.pkl'
    cfg_pkl = f'{workdir}/config.pkl'
    landscape_dir = f'{workdir}/landscape.{E}' if args.landscape_dir is None else args.landscape_dir
    outdir = f'{landscape_dir}/landscape_full' if args.outdir is None else args.outdir

    mask_mrc = f'{landscape_dir}/mask.mrc'
    pca_obj_pkl = f'{landscape_dir}/vol_pca_obj.pkl'
    assert os.path.exists(mask_mrc), f"{mask_mrc} missing. Did you run cryodrgn analyze_landscape?"
    assert os.path.exists(pca_obj_pkl), f"{pca_obj_pkl} missing. Did you run cryodrgn analyze_landscape?"
    
    log(f'Saving results to {outdir}')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    embeddings_pkl = f'{outdir}/vol_pca_sampled.pkl'
    z_sampled_pkl = f'{outdir}/z.sampled.pkl'
    if args.skip_vol:
        assert os.path.exists(embeddings_pkl), f"{embeddings_pkl} missing. Can't use --skip-vol"
        assert os.path.exists(z_sampled_pkl), f"{z_sampled_pkl} missing. Can't use --skip-vol"
        embeddings = utils.load_pkl(embeddings_pkl).astype(np.float32)
        z = utils.load_pkl(z_sampled_pkl)
    else:
        z, embeddings = generate_and_map_volumes(zfile, cfg_pkl, weights, mask_mrc, pca_obj_pkl,
                                                 landscape_dir, outdir, args)
        utils.save_pkl(embeddings, embeddings_pkl)

    # Train model
    embeddings_all = train_model(z, embeddings, outdir, zfile, args)
    utils.save_pkl(embeddings_all, f'{outdir}/vol_pca_all.pkl')

    # Copy viz notebook
    out_ipynb = f'{landscape_dir}/cryoDRGN_analyze_landscape.ipynb'
    if not os.path.exists(out_ipynb):
        log(f'Creating jupyter notebook...')
        ipynb = f'{cryodrgn._ROOT}/templates/cryoDRGN_analyze_landscape_template.ipynb'
        shutil.copyfile(ipynb, out_ipynb)
    else:
        log(f'{out_ipynb} already exists. Skipping')
    log(out_ipynb)
   
    log(f'Finished in {dt.now()-t1}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)

