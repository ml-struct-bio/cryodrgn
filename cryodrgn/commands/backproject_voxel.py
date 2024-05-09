"""Backproject cryo-EM images to reconstruct a volume as well as half-maps.

This script creates a reconstructed volume using the images in the given stack as well
as the given poses. Unless instructed otherwise, it will also produce volumes using
the images in each half of the dataset, as well as calculating an FSC curve between
these two half-map reconstructions.

Example usages
----------
$ cryodrgn backproject_voxel particles.128.mrcs --poses pose.pkl -o backproj.128.mrc
$ cryodrgn backproject_voxel particles.256.mrcs --poses pose.pkl
                             --ind good-particles.pkl -o backproj.256.mrc --lazy

"""
import argparse
import os
import time
import numpy as np
import torch
import matplotlib.pyplot as plt
import logging

from cryodrgn import ctf, dataset, fft, utils
from cryodrgn.mrc import MRCFile
from cryodrgn.lattice import Lattice
from cryodrgn.pose import PoseTracker
from cryodrgn.commands_utils.fsc import calculate_fsc

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
    parser.add_argument(
        "--no-half-maps",
        action="store_false",
        help="Don't produce half-maps and FSCs.",
        dest="half_maps",
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


def add_slice(volume, counts, ff_coord, ff, D, ctf_mul):
    d2 = int(D / 2)
    ff_coord = ff_coord.transpose(0, 1)
    xf, yf, zf = ff_coord.floor().long()
    xc, yc, zc = ff_coord.ceil().long()

    def add_for_corner(xi, yi, zi):
        dist = torch.stack([xi, yi, zi]).float() - ff_coord
        w = 1 - dist.pow(2).sum(0).pow(0.5)
        w[w < 0] = 0
        volume[(zi + d2, yi + d2, xi + d2)] += w * ff * ctf_mul
        counts[(zi + d2, yi + d2, xi + d2)] += w * ctf_mul**2

    add_for_corner(xf, yf, zf)
    add_for_corner(xc, yf, zf)
    add_for_corner(xf, yc, zf)
    add_for_corner(xf, yf, zc)
    add_for_corner(xc, yc, zf)
    add_for_corner(xf, yc, zc)
    add_for_corner(xc, yf, zc)
    add_for_corner(xc, yc, zc)


def regularize_volume(volume, counts, reg_weight):
    regularized_counts = counts + reg_weight * counts.mean()
    regularized_counts *= counts.mean() / regularized_counts.mean()
    reg_volume = volume / regularized_counts

    return fft.ihtn_center(reg_volume[0:-1, 0:-1, 0:-1].cpu())


def main(args):
    if not args.o.endswith(".mrc"):
        raise ValueError(f"Output file {args.o} does not end with .mrc!")

    t1 = time.time()
    logger.info(args)
    os.makedirs(os.path.dirname(args.o), exist_ok=True)

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
        assert (
            args.dose_per_tilt is not None
        ), "Argument --dose-per-tilt is required for backprojecting tilt series data"
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
        logger.info(f"Loading ctf params from {args.ctf}")
        ctf_params = ctf.load_ctf_for_training(D - 1, args.ctf)

        if args.ind is not None:
            ctf_params = ctf_params[args.ind]
        if args.tilt:
            ctf_params = np.concatenate(
                (ctf_params, data.ctfscalefactor.reshape(-1, 1)), axis=1  # type: ignore
            )
        ctf_params = torch.tensor(ctf_params, device=device)

    else:
        ctf_params = None

    Apix = float(ctf_params[0, 0]) if ctf_params is not None else 1.0
    voltage = float(ctf_params[0, 4]) if ctf_params is not None else None
    data.voltage = voltage
    mask = lattice.get_circular_mask(D // 2)
    iterator = range(min(args.first, Nimg)) if args.first else range(Nimg)

    volume_full = torch.zeros((D, D, D), device=device)
    counts_full = torch.zeros((D, D, D), device=device)
    volume_half1 = torch.zeros((D, D, D), device=device)
    counts_half1 = torch.zeros((D, D, D), device=device)
    volume_half2 = torch.zeros((D, D, D), device=device)
    counts_half2 = torch.zeros((D, D, D), device=device)

    for ii in iterator:
        if ii % 100 == 0:
            logger.info(f"fimage {ii}")

        r, t = posetracker.get_pose(ii)
        ff = data.get_tilt(ii) if args.tilt else data[ii]
        assert isinstance(ff, tuple)

        ff = ff[0].to(device)
        ff = ff.view(-1)[mask]
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
        add_slice(volume_full, counts_full, ff_coord, ff, D, ctf_mul)

        if args.half_maps:
            if ii % 2 == 0:
                add_slice(volume_half1, counts_half1, ff_coord, ff, D, ctf_mul)
            else:
                add_slice(volume_half2, counts_half2, ff_coord, ff, D, ctf_mul)

    td = time.time() - t1
    logger.info(
        f"Backprojected {len(iterator)} images "
        f"in {td:.2f}s ({(td / Nimg):4f}s per image)"
    )

    counts_full[counts_full == 0] = 1
    counts_half1[counts_half1 == 0] = 1
    counts_half2[counts_half2 == 0] = 1

    if args.output_sumcount:
        MRCFile.write(args.o + ".sums", volume_full.cpu().numpy(), Apix=Apix)
        MRCFile.write(args.o + ".counts", counts_full.cpu().numpy(), Apix=Apix)

    volume_full = regularize_volume(volume_full, counts_full, args.reg_weight)
    out_path = os.path.splitext(args.o)[0]
    MRCFile.write(args.o, np.array(volume_full).astype("float32"), Apix=Apix)

    if args.half_maps:
        volume_half1 = regularize_volume(volume_half1, counts_half1, args.reg_weight)
        volume_half2 = regularize_volume(volume_half2, counts_half2, args.reg_weight)
        fsc_vals = calculate_fsc(volume_half1, volume_half2)

        fsc_vals.to_csv("_".join([out_path, "fsc-vals.txt"]), sep=" ", header=False)
        plt.plot(fsc_vals.pixres, fsc_vals.fsc)
        plt.ylim((0, 1))
        plt.savefig("_".join([out_path, "fsc-plot.png"]), bbox_inches="tight")

        if ((fsc_vals.fsc >= 0.5) & (fsc_vals.pixres > 0)).any():
            fsc_res = fsc_vals.pixres[fsc_vals.fsc >= 0.5].max() ** -1.0 * Apix
            logger.info(f"res @ FSC=0.5: {fsc_res:.4f}")
        else:
            logger.warning("res @ FSC=0.5: N/A")
        if ((fsc_vals.fsc >= 0.143) & (fsc_vals.pixres > 0)).any():
            fsc_res = fsc_vals.pixres[fsc_vals.fsc >= 0.143].max() ** -1.0 * Apix
            logger.info(f"res @ FSC=0.143: {fsc_res:.4f}")
        else:
            logger.warning("res @ FSC=0.143: N/A")

        # save the half-map reconstructions to file
        half_fl1 = "_".join([out_path, "half-map1.mrc"])
        half_fl2 = "_".join([out_path, "half-map2.mrc"])
        MRCFile.write(half_fl1, np.array(volume_half1).astype("float32"), Apix=Apix)
        MRCFile.write(half_fl2, np.array(volume_half2).astype("float32"), Apix=Apix)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
