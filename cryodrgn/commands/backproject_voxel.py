"""
Backproject cryo-EM images
"""

import argparse
import os
import time
import numpy as np
import torch
import logging
from cryodrgn import ctf, dataset, fft, utils
from cryodrgn.mrc import MRCFile
from cryodrgn.lattice import Lattice
from cryodrgn.pose import PoseTracker

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument(
        "particles",
        type=os.path.abspath,
        help="Input particles (.mrcs, .star, .cs, or .txt)",
    )
    parser.add_argument(
        "--poses", type=os.path.abspath, required=True, help="Image poses (.pkl)"
    )
    parser.add_argument(
        "--ctf",
        metavar="pkl",
        type=os.path.abspath,
        help="CTF parameters (.pkl) for phase flipping images",
    )
    parser.add_argument(
        "-o", type=os.path.abspath, required=True, help="Output .mrc file"
    )
    parser.add_argument("--ctf-alg", type=str, choices=("flip", "mul"), default="mul")
    parser.add_argument(
        "--reg-weight",
        type=float,
        default=1.0,
        help="Add this value times the mean weight to the weight map to regularize the volume, reducing noise."
        "Alternatively, you can set --output-sumcount, and then use `cryodrgn_utils regularize_backproject` on the"
        ".sums and .counts files to try different regularization constants post hoc.",
    )
    parser.add_argument(
        "--output-sumcount",
        action="store_true",
        help="Output voxel sums and counts so that different regularization weights can be applied post hoc, with "
        "`cryodrgn_utils regularize_backproject`.",
    )

    group = parser.add_argument_group("Dataset loading options")
    group.add_argument(
        "--uninvert-data",
        dest="invert_data",
        action="store_false",
        help="Do not invert data sign",
    )
    group.add_argument(
        "--datadir",
        type=os.path.abspath,
        help="Path prefix to particle stack if loading relative paths from a .star or .cs file",
    )
    group.add_argument(
        "--lazy",
        action="store_true",
        help="Lazy loading if full dataset is too large to fit in memory",
    )
    group.add_argument(
        "--ind",
        type=os.path.abspath,
        metavar="PKL",
        help="Filter particles by these indices",
    )
    group.add_argument(
        "--first",
        type=int,
        default=None,
        help="Backproject the first N images (default: all images)",
    )
    group = parser.add_argument_group("Tilt series options")
    group.add_argument(
        "--tilt",
        action="store_true",
        help="Flag to treat data as a tilt series from cryo-ET",
    )
    group.add_argument(
        "--ntilts",
        type=int,
        default=10,
        help="Number of tilts per particle to backproject (default: %(default)s)",
    )
    group.add_argument(
        "-d",
        "--dose-per-tilt",
        type=float,
        help="Expected dose per tilt (electrons/A^2 per tilt) (default: %(default)s)",
    )
    group.add_argument(
        "-a",
        "--angle-per-tilt",
        type=float,
        default=3,
        help="Tilt angle increment per tilt in degrees (default: %(default)s)",
    )

    return parser


def add_slice(V, counts, ff_coord, ff, D, ctf_mul):
    d2 = int(D / 2)
    ff_coord = ff_coord.transpose(0, 1)
    xf, yf, zf = ff_coord.floor().long()
    xc, yc, zc = ff_coord.ceil().long()

    def add_for_corner(xi, yi, zi):
        dist = torch.stack([xi, yi, zi]).float() - ff_coord
        w = 1 - dist.pow(2).sum(0).pow(0.5)
        w[w < 0] = 0
        V[(zi + d2, yi + d2, xi + d2)] += w * ff * ctf_mul
        counts[(zi + d2, yi + d2, xi + d2)] += w * ctf_mul**2

    add_for_corner(xf, yf, zf)
    add_for_corner(xc, yf, zf)
    add_for_corner(xf, yc, zf)
    add_for_corner(xf, yf, zc)
    add_for_corner(xc, yc, zf)
    add_for_corner(xf, yc, zc)
    add_for_corner(xc, yf, zc)
    add_for_corner(xc, yc, zc)
    return V, counts


def main(args):
    assert args.o.endswith(".mrc")

    t1 = time.time()
    logger.info(args)
    if not os.path.exists(os.path.dirname(args.o)):
        os.makedirs(os.path.dirname(args.o))

    # set the device
    use_cuda = torch.cuda.is_available()
    device = torch.device("cuda" if use_cuda else "cpu")
    logger.info("Use cuda {}".format(use_cuda))
    if not use_cuda:
        logger.warning("WARNING: No GPUs detected")

    # load the particles
    if args.ind is not None:
        if args.tilt:
            particle_ind = utils.load_pkl(args.ind).astype(int)
            pt, tp = dataset.TiltSeriesData.parse_particle_tilt(args.particles)
            tilt_ind = dataset.TiltSeriesData.particles_to_tilts(pt, particle_ind)
            args.ind = tilt_ind
        else:
            args.ind = utils.load_pkl(args.ind).astype(int)

    if args.tilt:
        data = dataset.TiltSeriesData(
            args.particles,
            args.ntilts,
            norm=(0, 1),
            invert_data=args.invert_data,
            datadir=args.datadir,
            ind=args.ind,
            lazy=args.lazy,
            dose_per_tilt=args.dose_per_tilt,
            angle_per_tilt=args.angle_per_tilt,
            device=device,
        )
    else:
        data = dataset.ImageDataset(
            mrcfile=args.particles,
            norm=(0, 1),
            invert_data=args.invert_data,
            datadir=args.datadir,
            ind=args.ind,
            lazy=args.lazy,
        )

    D = data.D
    Nimg = data.N

    lattice = Lattice(D, extent=D // 2, device=device)

    posetracker = PoseTracker.load(args.poses, Nimg, D, None, args.ind, device=device)

    if args.ctf is not None:
        logger.info("Loading ctf params from {}".format(args.ctf))
        ctf_params = ctf.load_ctf_for_training(D - 1, args.ctf)
        if args.ind is not None:
            ctf_params = ctf_params[args.ind]
        if args.tilt:
            assert (
                args.dose_per_tilt is not None
            ), "Argument --dose-per-tilt is required for backprojecting tilt series data"
            ctf_params = np.concatenate(
                (ctf_params, data.ctfscalefactor.reshape(-1, 1)), axis=1  # type: ignore
            )
        ctf_params = torch.tensor(ctf_params, device=device)
    else:
        ctf_params = None
    Apix = float(ctf_params[0, 0]) if ctf_params is not None else 1.0
    voltage = float(ctf_params[0, 4]) if ctf_params is not None else None
    data.voltage = voltage

    V = torch.zeros((D, D, D), device=device)
    counts = torch.zeros((D, D, D), device=device)

    mask = lattice.get_circular_mask(D // 2)

    if args.first:
        args.first = min(args.first, Nimg)
        iterator = range(args.first)
    else:
        iterator = range(Nimg)

    for ii in iterator:
        if ii % 100 == 0:
            logger.info("image {}".format(ii))
        r, t = posetracker.get_pose(ii)
        ff = data.get_tilt(ii) if args.tilt else data[ii]
        assert isinstance(ff, tuple)

        ff = ff[0].to(device)
        ff = ff.view(-1)[mask]
        c = None
        ctf_mul = 1
        if ctf_params is not None:
            freqs = lattice.freqs2d / ctf_params[ii, 0]
            c = ctf.compute_ctf(freqs, *ctf_params[ii, 1:]).view(-1)[mask]
            if args.ctf_alg == "flip":
                ff *= c.sign()
            else:
                ctf_mul = c
        if t is not None:
            ff = lattice.translate_ht(ff.view(1, -1), t.view(1, 1, 2), mask).view(-1)
        if args.tilt:
            tilt_idxs = torch.tensor([ii]).to(device)
            dose_filters = data.get_dose_filters(tilt_idxs, lattice, ctf_params[ii, 0])[
                0
            ]
            ctf_mul *= dose_filters[mask]

        ff_coord = lattice.coords[mask] @ r
        add_slice(V, counts, ff_coord, ff, D, ctf_mul)

    td = time.time() - t1
    logger.info(
        "Backprojected {} images in {}s ({}s per image)".format(
            len(iterator), td, td / Nimg
        )
    )
    counts[counts == 0] = 1

    if args.output_sumcount:
        MRCFile.write(args.o + ".sums", V.cpu().numpy(), Apix=Apix)
        MRCFile.write(args.o + ".counts", counts.cpu().numpy(), Apix=Apix)

    regularized_counts = counts + args.reg_weight * counts.mean()
    regularized_counts *= counts.mean() / regularized_counts.mean()
    V /= regularized_counts
    V = fft.ihtn_center(V[0:-1, 0:-1, 0:-1].cpu())
    MRCFile.write(args.o, np.array(V).astype("float32"), Apix=Apix)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    main(add_args(parser).parse_args())
