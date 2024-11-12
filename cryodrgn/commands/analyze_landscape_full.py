"""Transform a cryoDRGN latent space to better capture differences between volumes.

Example usage
-------------
$ cryodrgn analyze_landscape_full 003_abinit-het/ 49

# Sample fewer volumes from this dataset's particles according to their position in the
# latent space; use a larger box size for these volumes instead of downsampling to 128
$ cryodrgn analyze_landscape_full 005_train-vae/ 39 -N 4000 -d 256

"""
import argparse
import os
import pprint
import shutil
from datetime import datetime as dt
import logging
import nbformat

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.model_selection import train_test_split
from torch.utils.data.dataloader import default_collate
from torch.utils.data import Dataset, DataLoader

import cryodrgn
from cryodrgn import config, utils
from cryodrgn.models.neural_nets import ResidLinearMLP
from cryodrgn.models.utils import load_model
from cryodrgn.source import ImageSource

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments for use with `cryodrgn analyze_landscape_full`."""

    parser.add_argument(
        "workdir", type=os.path.abspath, help="Directory with cryoDRGN results"
    )
    parser.add_argument(
        "epoch",
        type=int,
        help="Epoch number N to analyze (0-based indexing, "
        "corresponding to z.N.pkl and weights.N.pkl)",
    )
    parser.add_argument("--device", type=int, help="Optionally specify CUDA device")
    parser.add_argument(
        "--landscape-dir",
        type=os.path.abspath,
        help="Landscape analysis directory (default: [workdir]/landscape.[epoch])",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=os.path.abspath,
        help="Output directory (default: [workdir]/landscape.[epoch]/landscape_full)",
    )
    parser.add_argument(
        "--seed", default=0, type=int, help="Random seed (default: %(default)s)"
    )

    group = parser.add_argument_group("Volume generation arguments")
    group.add_argument(
        "-N",
        "--training-volumes",
        type=int,
        default=10000,
        help="Number of training volumes to generate (default: %(default)s)",
    )
    group.add_argument("--flip", action="store_true", help="Flip handedness")
    group.add_argument(
        "-d",
        "--downsample",
        type=int,
        default=128,
        help="Downsample volumes to this box size (pixels) (default: %(default)s)",
    )
    group.add_argument(
        "--skip-vol", action="store_true", help="Skip generation of volumes"
    )

    group = parser.add_argument_group("Volume mapping arguments")
    group.add_argument(
        "--batch-size",
        type=int,
        default=64,
        metavar="N",
        help="input batch size for training (default: 64)",
    )
    group.add_argument(
        "--test-batch-size",
        type=int,
        default=1000,
        metavar="N",
        help="input batch size for testing (default: 1000)",
    )
    group.add_argument(
        "--epochs",
        type=int,
        default=50,
        metavar="N",
        help="number of epochs to train (default: 14)",
    )
    group.add_argument(
        "--lr",
        type=float,
        default=1e-3,
        metavar="LR",
        help="learning rate (default: 1e-3)",
    )
    group.add_argument(
        "--dim",
        type=int,
        default=512,
        metavar="N",
        help="MLP hidden layer dimension (default: 512)",
    )
    group.add_argument(
        "--layers",
        type=int,
        default=3,
        metavar="N",
        help="MLP number of hidden layers (default: 3)",
    )


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
            print(
                "Train Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}".format(
                    epoch,
                    batch_idx * len(data),
                    len(train_loader.dataset),
                    100.0 * batch_idx / len(train_loader),
                    loss.item(),
                )
            )


def test(model, device, test_loader):
    model.eval()
    test_loss = 0
    with torch.no_grad():
        for data, target in test_loader:
            data, target = data.to(device), target.to(device)
            output = model(data)
            test_loss += (
                len(data) * F.mse_loss(output, target).item()
            )  # sum up batch loss

    test_loss /= len(test_loader.dataset)

    print("\nTest set: Average loss: {:.4f}\n".format(test_loss))


class MyDataset(Dataset):
    def __init__(self, x, y):
        assert len(x) == len(y)
        self.x = x
        self.y = y

    def __len__(self):
        return len(self.x)

    def __getitem__(self, idx):
        return self.x[idx], self.y[idx]


def generate_and_map_volumes(
    zfile, cfg, weights, mask_mrc, pca_obj_pkl, landscape_dir, outdir, args
):
    # Sample z
    logger.info(f"Sampling {args.training_volumes} particles from {zfile}")
    np.random.seed(args.seed)
    z_all = utils.load_pkl(zfile)
    ind = np.array(
        sorted(np.random.choice(len(z_all), args.training_volumes, replace=False))
    )  # type: ignore
    z_sample = z_all[ind]
    utils.save_pkl(z_sample, f"{outdir}/z.sampled.pkl")
    utils.save_pkl(ind, f"{outdir}/ind.sampled.pkl")
    logger.info(f"Saved {outdir}/z.sampled.pkl")

    # Set the device
    # if torch.cuda.is_available():
    #     torch.set_default_tensor_type(torch.cuda.FloatTensor)  # type: ignore

    cfg = config.update_config_v1(cfg)
    logger.info("Loaded configuration:")
    pprint.pprint(cfg)

    D = cfg["lattice_args"]["D"]  # image size + 1
    norm = [float(x) for x in cfg["dataset_args"]["norm"]]

    # Load landscape analysis inputs
    mask = np.array(ImageSource.from_file(mask_mrc).images().cpu())
    assert isinstance(mask, np.ndarray)
    mask = mask.astype(bool)
    if args.downsample:
        assert mask.shape == (args.downsample,) * 3
    else:
        assert mask.shape == (D - 1, D - 1, D - 1)
    logger.info(f"{mask.sum()} voxels in the mask")

    pca = utils.load_pkl(pca_obj_pkl)

    # Load model weights
    logger.info("Loading weights from {}".format(weights))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, lattice = load_model(cfg, weights, device)
    model.eval()

    # Set z
    z = z_sample.astype(np.float32)

    # Generate volumes
    logger.info(f"Generating {len(z)} volume embeddings")
    t1 = dt.now()
    embeddings = []
    for i, zz in enumerate(z):
        if i % 100 == 0:
            logger.info(i)

        if args.downsample:
            extent = lattice.extent * (args.downsample / (D - 1))
            vol = model.eval_volume(
                lattice.get_downsample_coords(args.downsample + 1),
                args.downsample + 1,
                extent,
                norm,
                zz,
            )
        else:
            vol = model.eval_volume(lattice.coords, lattice.D, lattice.extent, norm, zz)
        if args.flip:
            vol = vol.flip([0])

        embeddings.append(
            pca.transform(vol.cpu()[torch.tensor(mask).bool()].reshape(1, -1))
        )

    embeddings = np.array(embeddings).reshape(len(z), -1).astype(np.float32)
    td = dt.now() - t1
    logger.info(f"Finished generating {args.training_volumes} volumes in {td}")

    return z, embeddings


def train_model(x, y, outdir, zfile, args):
    use_cuda = torch.cuda.is_available()
    torch.manual_seed(args.seed)
    device = torch.device("cuda" if use_cuda else "cpu")

    train_kwargs = {"batch_size": args.batch_size}
    test_kwargs = {"batch_size": args.test_batch_size}
    if use_cuda:
        cuda_kwargs = {
            "num_workers": 0,
            "shuffle": True,
            "collate_fn": lambda x: tuple(x_.to(device) for x_ in default_collate(x)),
        }
        train_kwargs.update(cuda_kwargs)
        test_kwargs.update(cuda_kwargs)

    # Load dataset
    x_train, x_test, y_train, y_test = train_test_split(
        x, y, test_size=0.25, random_state=args.seed
    )

    train_dataset = MyDataset(x_train, y_train)
    test_dataset = MyDataset(x_test, y_test)

    train_loader = DataLoader(train_dataset, **train_kwargs)
    test_loader = DataLoader(test_dataset, **test_kwargs)

    model = ResidLinearMLP(x.shape[1], args.layers, args.dim, y.shape[1], nn.ReLU).to(
        device
    )
    logger.info(model)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    # Train
    for epoch in range(1, args.epochs + 1):
        train(args, model, device, train_loader, optimizer, epoch)
        test(model, device, test_loader)

    # Evaluate
    model.eval()
    yhat_all = []
    eval_dataset = utils.load_pkl(zfile).astype(np.float32)
    with torch.no_grad():
        for x in np.array_split(eval_dataset, args.test_batch_size):
            x = torch.tensor(x, device=device)
            yhat = model(x)
            yhat_all.append(yhat.detach().cpu().numpy())

    yhat_all = np.concatenate(yhat_all)
    torch.save(model.state_dict(), f"{outdir}/model.pt")

    return yhat_all


def main(args: argparse.Namespace) -> None:
    t1 = dt.now()
    logger.info(args)

    E = args.epoch
    workdir = args.workdir
    zfile = f"{workdir}/z.{E}.pkl"
    weights = f"{workdir}/weights.{E}.pkl"
    cfg = (
        f"{workdir}/config.yaml"
        if os.path.exists(f"{workdir}/config.yaml")
        else f"{workdir}/config.pkl"
    )
    landscape_dir = (
        f"{workdir}/landscape.{E}" if args.landscape_dir is None else args.landscape_dir
    )
    outdir = f"{landscape_dir}/landscape_full" if args.outdir is None else args.outdir

    mask_mrc = f"{landscape_dir}/mask.mrc"
    pca_obj_pkl = f"{landscape_dir}/vol_pca_obj.pkl"
    assert os.path.exists(
        mask_mrc
    ), f"{mask_mrc} missing. Did you run cryodrgn analyze_landscape?"
    assert os.path.exists(
        pca_obj_pkl
    ), f"{pca_obj_pkl} missing. Did you run cryodrgn analyze_landscape?"

    kmeans_folder = [
        p for p in os.listdir(landscape_dir) if p.startswith("clustering_L2_")
    ]
    if len(kmeans_folder) == 0:
        raise RuntimeError(
            "No clustering folders `clustering_L2_` found. "
            "Did you run cryodrgn analyze_landscape?"
        )
    kmeans_folder = kmeans_folder[0]
    link_method, kmeans_K = kmeans_folder.split("_")[-2:]

    logger.info(f"Saving results to {outdir}")
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    embeddings_pkl = f"{outdir}/vol_pca_sampled.pkl"
    z_sampled_pkl = f"{outdir}/z.sampled.pkl"
    if args.skip_vol:
        assert os.path.exists(
            embeddings_pkl
        ), f"{embeddings_pkl} missing. Can't use --skip-vol"
        assert os.path.exists(
            z_sampled_pkl
        ), f"{z_sampled_pkl} missing. Can't use --skip-vol"
        embeddings = utils.load_pkl(embeddings_pkl).astype(np.float32)
        z = utils.load_pkl(z_sampled_pkl)
    else:
        z, embeddings = generate_and_map_volumes(
            zfile, cfg, weights, mask_mrc, pca_obj_pkl, landscape_dir, outdir, args
        )
        utils.save_pkl(embeddings, embeddings_pkl)

    # Train model
    embeddings_all = train_model(z, embeddings, outdir, zfile, args)
    utils.save_pkl(embeddings_all, f"{outdir}/vol_pca_all.pkl")

    # Copy viz notebook
    out_ipynb = os.path.join(landscape_dir, "cryoDRGN_analyze_landscape.ipynb")
    if not os.path.exists(out_ipynb):
        logger.info("Creating jupyter notebook...")
        ipynb = os.path.join(
            cryodrgn._ROOT, "templates", "cryoDRGN_analyze_landscape_template.ipynb"
        )
        shutil.copyfile(ipynb, out_ipynb)
    else:
        logger.info(f"{out_ipynb} already exists. Skipping")

    # Lazily look at the beginning of the notebook for the epoch number to update
    with open(out_ipynb, "r") as f:
        filter_ntbook = nbformat.read(f, as_version=nbformat.NO_CONVERT)

    for cell in filter_ntbook["cells"]:
        cell["source"] = cell["source"].replace("EPOCH = None", f"EPOCH = {args.epoch}")
        cell["source"] = cell["source"].replace(
            "WORKDIR = None", f'WORKDIR = "{args.workdir}"'
        )
        cell["source"] = cell["source"].replace(
            "K = None", f"K = {args.training_volumes}"
        )
        cell["source"] = cell["source"].replace("M = None", f"M = {kmeans_K}")
        cell["source"] = cell["source"].replace(
            "linkage = None", f'linkage = "{link_method}"'
        )

    with open(out_ipynb, "w") as f:
        nbformat.write(filter_ntbook, f)

    logger.info(out_ipynb)
    logger.info(f"Finished in {dt.now()-t1}")
