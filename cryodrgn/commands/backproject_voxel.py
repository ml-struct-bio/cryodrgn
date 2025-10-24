"""Voxel-based backprojection to reconstruct a volume as well as half-maps.

This command performs volume reconstruction on the given image stack using voxel-based
backprojection applied to as well as the given poses. Unless instructed otherwise,
it will also produce volumes using the images in each half of the dataset, along with
calculating an FSC curve between these two half-map reconstructions.
All outputs will be placed in the folder specified.

Example usage
-------------
$ cryodrgn backproject_voxel particles.128.mrcs \
                             --ctf ctf.pkl --poses pose.pkl -o backproject-results/

# Use `--lazy` for large datasets to reduce memory usage and avoid OOM errors
$ cryodrgn backproject_voxel particles.256.mrcs --ctf ctf.pkl --poses pose.pkl \
                             --ind good_particles.pkl -o backproject-256/ --lazy

# `--first` is also a good tool for doing a quick initial reconstruction with only the
# first <x> images in the given stack
# Here we also avoid calculating FSCs between the half-maps at the end
$ cryodrgn backproject_voxel particles.196.mrcs --ctf ctf.pkl --poses pose.pkl \
                             -o backproject-196/ --lazy --first=10000

# `--tilt` is required for subtomogram datasets; you can further control how many tilts
# are used per particle using `--ntilts`
# `--datadir` is generally required when using .star or .cs particle inputs
$ cryodrgn backproject_voxel particles_from_M.star --datadir subtilts/128/ \
                             --ctf ctf.pkl --poses pose.pkl \
                             -o backproject-tilt/ --lazy --tilt --ntilts 5

"""
import argparse
import os
import time
import numpy as np
import torch
import logging
from cryodrgn import ctf, dataset, fft, utils
from cryodrgn.lattice import Lattice
from cryodrgn.pose import PoseTracker
from cryodrgn.commands_utils.fsc import calculate_cryosparc_fscs
from cryodrgn.source import write_mrc

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments for use with `cryodrgn backproject_voxel`."""

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
        "--outdir",
        "-o",
        type=os.path.abspath,
        required=True,
        help="New or existing folder in which outputs will be " "placed",
    )
    parser.add_argument(
        "--no-half-maps",
        action="store_false",
        help="Don't produce half-maps and FSCs.",
        dest="half_maps",
    )
    parser.add_argument(
        "--no-fsc-vals",
        action="store_false",
        help="Don't calculate FSCs, but still produce half-maps.",
        dest="fsc_vals",
    )
    parser.add_argument(
        "-b",
        "--batch-size",
        type=int,
        default=1000,
        help="Number of images to process in each batch (default: %(default)s)",
    )

    parser.add_argument("--ctf-alg", type=str, choices=("flip", "mul"), default="mul")
    parser.add_argument(
        "--reg-weight",
        type=float,
        default=1.0,
        help="Add this value times the mean weight to the weight map to regularize the"
        "volume, reducing noise.\nAlternatively, you can set --output-sumcount, and "
        "then use `cryodrgn_utils regularize_backproject` on the"
        ".sums and .counts files to try different regularization constants post hoc.\n"
        "(default: %(default)s)",
    )
    parser.add_argument(
        "--output-sumcount",
        action="store_true",
        help="Output voxel sums and counts so that different regularization weights "
        "can be applied post hoc, with `cryodrgn_utils regularize_backproject`.",
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
        help="Path prefix to particle stack files if loading "
        "relative stack paths from a .star or .cs file",
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
        help="Filter particles by these indices before starting backprojection",
    )
    group.add_argument(
        "--first",
        type=int,
        help="Backproject using only the first N images (default: all images)",
    )
    parser.add_argument(
        "--log-interval",
        type=str,
        default="5000",
        help="Logging interval in N_IMGS (default: %(default)s)",
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
        "--force-ntilts",
        action="store_true",
        dest="force_ntilts",
        help="Automatically keep only particles with ≥ --ntilts tilts",
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
    res = volume.shape[0]
    use_coord = ff_coord.transpose(0, 2).clip(-d2, d2)
    xf, yf, zf = use_coord.floor().long()
    xc, yc, zc = use_coord.ceil().long()

    def add_for_corner(xi, yi, zi):
        dist = torch.stack([xi, yi, zi]).float() - use_coord
        wvals = 1 - dist.pow(2).sum(0).pow(0.5).T
        wvals[wvals < 0] = 0
        vol_add = wvals * ff * ctf_mul
        cnt_add = wvals * ctf_mul**2
        flat_idx = (
            (zi.T.flatten() + d2) * res**2
            + (yi.T.flatten() + d2) * res
            + (xi.T.flatten() + d2)
        )
        volume.put_(flat_idx, vol_add.flatten(), accumulate=True)
        counts.put_(flat_idx, cnt_add.flatten(), accumulate=True)

    add_for_corner(xf, yf, zf)
    add_for_corner(xc, yf, zf)
    add_for_corner(xf, yc, zf)
    add_for_corner(xf, yf, zc)
    add_for_corner(xc, yc, zf)
    add_for_corner(xf, yc, zc)
    add_for_corner(xc, yf, zc)
    add_for_corner(xc, yc, zc)


def regularize_volume(
    volume: torch.Tensor, counts: torch.Tensor, reg_weight: float
) -> torch.Tensor:
    regularized_counts = counts + reg_weight * counts.mean()
    regularized_counts *= counts.mean() / regularized_counts.mean()
    reg_volume = volume / regularized_counts

    return fft.ihtn_center(reg_volume[0:-1, 0:-1, 0:-1].cpu())


def main(args: argparse.Namespace) -> None:
    t1 = time.time()
    logger.info(args)
    os.makedirs(args.outdir, exist_ok=True)

    # Set the device for doing computation on (GPU if available)
    use_cuda = torch.cuda.is_available()
    device = torch.device("cuda" if use_cuda else "cpu")
    logger.info("Use cuda {}".format(use_cuda))
    if not use_cuda:
        logger.warning("WARNING: No GPUs detected")

    # Generate indices if forcing ntilts and save them to the given file
    if args.force_ntilts:
        if not args.tilt:
            raise ValueError("--force-ntilts is only supported for tilt series data!")
        if args.ind is not None:
            raise ValueError("--force-ntilts overrides --ind!")
        logger.info(
            f"Filtering to particles with ≥ {args.ntilts} tilts (–force-ntilts)"
        )

        pt, _ = dataset.TiltSeriesData.parse_particle_tilt(args.particles)
        counts = [len(tilt_list) for tilt_list in pt]
        valid_particles = np.where(np.array(counts) >= args.ntilts)[0]
        if valid_particles.size == 0:
            raise ValueError(f"No particles have at least {args.ntilts} tilts.")
        idx_file = os.path.join(args.outdir, f"indices_force_ntilts_{args.ntilts}.pkl")
        utils.save_pkl(valid_particles.astype(int), idx_file)
        logger.info(f"→ saved {valid_particles.size} particle IDs to {idx_file}")
        args.ind = idx_file

    # Load the filtering indices and the image dataset to use for backprojection
    if args.ind is not None:
        if args.tilt:
            particle_ind = utils.load_pkl(args.ind).astype(int)
            pt, tp = dataset.TiltSeriesData.parse_particle_tilt(args.particles)
            tilt_ind = dataset.TiltSeriesData.particles_to_tilts(pt, particle_ind)
            indices = tilt_ind
        else:
            indices = utils.load_pkl(args.ind).astype(int)
    else:
        indices = None

    if args.tilt:
        assert (
            args.dose_per_tilt is not None
        ), "Argument --dose-per-tilt is required for backprojecting tilt series data"
        data = dataset.TiltSeriesData(
            args.particles,
            ntilts=None,
            norm=(0, 1),
            invert_data=args.invert_data,
            datadir=args.datadir,
            ind=indices,
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
            ind=indices,
            lazy=args.lazy,
        )

    D = data.D
    lattice = Lattice(D, extent=D // 2, device=device)
    posetracker = PoseTracker.load(args.poses, data.N, D, None, indices, device=device)

    if args.ctf is not None:
        logger.info(f"Loading ctf params from {args.ctf}")
        ctf_params = ctf.load_ctf_for_training(D - 1, args.ctf)

        if indices is not None:
            ctf_params = ctf_params[indices]
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
    lattice_mask = lattice.get_circular_mask(D // 2)
    mask_size = lattice_mask.sum().item()

    # Initialize tensors that will store backprojection results
    volume_full = torch.zeros((D, D, D), device=device)
    counts_full = torch.zeros((D, D, D), device=device)
    volume_half1 = torch.zeros((D, D, D), device=device)
    counts_half1 = torch.zeros((D, D, D), device=device)
    volume_half2 = torch.zeros((D, D, D), device=device)
    counts_half2 = torch.zeros((D, D, D), device=device)

    # Determine which images to backproject if only using the first N images from the
    # stack or the best N tilts per particle according to CTF scalefactor
    img_iterator = range(data.N)
    if args.tilt:
        use_tilts = set(range(args.ntilts))
        img_iterator = [
            ii for ii in img_iterator if int(data.tilt_numbers[ii].item()) in use_tilts
        ]
    else:
        img_iterator = list(img_iterator)

    if args.first:
        img_iterator = img_iterator[: args.first]

    # Figure out how often we are going to report progress w.r.t. images processed
    Nimg = len(img_iterator)
    Npart = 0 if args.tilt else Nimg
    if args.log_interval == "auto":
        log_interval = max(round((Nimg // 100), -2), 100)
    elif args.log_interval.isnumeric():
        log_interval = int(args.log_interval)
    else:
        raise ValueError(
            f"Unrecognized argument for --log-interval: `{args.log_interval}`\n"
            f"Given value must be an integer or the label 'auto'!"
        )

    # Training loop over batches of images from the dataset
    img_count = 0
    half_odd = False
    while img_iterator:
        img_idxs = img_iterator[: args.batch_size]
        img_iterator = img_iterator[args.batch_size :]
        B = len(img_idxs)
        img_count += B
        if args.tilt:
            Npart += len(set([data.get_tilt_particle(img_idx) for img_idx in img_idxs]))

        if (img_count % log_interval) < B:
            logger.info(f"fimage {img_count:,d} — {(img_count / Nimg * 100):.1f}% done")

        ff = data.get_tilt(img_idxs) if args.tilt else data[img_idxs]
        assert isinstance(ff, tuple)
        ff = ff[0].to(device)
        ff = ff.view(-1, D**2).to(device)[:, lattice_mask]
        rot, trans = posetracker.get_pose(img_idxs)
        rot = rot.to(device)
        if trans is not None:
            trans = trans.to(device)

        if ctf_params is not None:
            freqs = lattice.freqs2d.unsqueeze(0).expand(
                B, *lattice.freqs2d.shape
            ) / ctf_params[img_idxs, 0].view(B, 1, 1)
            c = ctf.compute_ctf(
                freqs, *torch.split(ctf_params[img_idxs, 1:], 1, 1)
            ).view(B, -1)[:, lattice_mask]

            if args.ctf_alg == "flip":
                ff *= c.sign()
                ctf_mul = torch.tensor(np.tile(1, (B, mask_size)), device=device)
            else:
                ctf_mul = c

        else:
            ctf_mul = torch.tensor(np.tile(1, (B, mask_size)), device=device)

        if trans is not None:
            ff = lattice.translate_ht(ff, trans.unsqueeze(1), lattice_mask).view(B, -1)

        if args.tilt:
            for i in range(B):
                dose_filters = data.get_dose_filters([img_idxs[i]], lattice, Apix)[0]
                ctf_mul[i] *= dose_filters[lattice_mask]

        # Accumulate backprojection results for this batch to the variables storing the
        # final results for the full dataset as well as the half-maps if applicable
        ff_coord = lattice.coords[lattice_mask] @ rot
        add_slice(volume_full, counts_full, ff_coord, ff, D, ctf_mul)
        if args.half_maps:
            even_idx = 1 if half_odd else 0
            odd_idx = 0 if half_odd else 1
            if not half_odd or B > 1:
                add_slice(
                    volume_half1,
                    counts_half1,
                    ff_coord[even_idx::2],
                    ff[even_idx::2],
                    D,
                    ctf_mul[even_idx::2],
                )
            if half_odd or B > 1:
                add_slice(
                    volume_half2,
                    counts_half2,
                    ff_coord[odd_idx::2],
                    ff[odd_idx::2],
                    D,
                    ctf_mul[odd_idx::2],
                )

            # When the number of images in the batch is odd, the first image in the
            # next batch will accumulate to the other half-map and vice versa
            if B % 2 == 1:
                half_odd = not half_odd

    td = time.time() - t1
    img_avg = td / Nimg
    if args.tilt:
        logger.info(
            f"Backprojected {Nimg:,d} tilts from {Npart:,d} particles in {td:.2f}s "
            f"({img_avg:4f}s per tilt image)"
        )
    else:
        logger.info(
            f"Backprojected {Nimg:,d} images in {td:.2f}s ({img_avg:4f}s per image)"
        )

    counts_full[counts_full == 0] = 1
    counts_half1[counts_half1 == 0] = 1
    counts_half2[counts_half2 == 0] = 1
    if args.output_sumcount:
        sums_fl = os.path.join(args.outdir, "backproject.sums")
        counts_fl = os.path.join(args.outdir, "backproject.counts")
        write_mrc(sums_fl, volume_full.cpu().numpy(), Apix=Apix)
        write_mrc(counts_fl, counts_full.cpu().numpy(), Apix=Apix)

    volume_full = regularize_volume(volume_full, counts_full, args.reg_weight)
    vol_fl = os.path.join(args.outdir, "backproject.mrc")
    write_mrc(vol_fl, np.array(volume_full).astype("float32"), Apix=Apix)

    # Create the half-maps, calculate the FSC curve between them, and save both to file
    if args.half_maps:
        volume_half1 = regularize_volume(volume_half1, counts_half1, args.reg_weight)
        volume_half2 = regularize_volume(volume_half2, counts_half2, args.reg_weight)
        half_fl1 = os.path.join(args.outdir, "half_map_a.mrc")
        half_fl2 = os.path.join(args.outdir, "half_map_b.mrc")
        write_mrc(half_fl1, np.array(volume_half1).astype("float32"), Apix=Apix)
        write_mrc(half_fl2, np.array(volume_half2).astype("float32"), Apix=Apix)

        if args.fsc_vals:
            out_file = os.path.join(args.outdir, "fsc-vals.txt")
            plot_file = os.path.join(args.outdir, "fsc-plot.png")
            _ = calculate_cryosparc_fscs(
                volume_full,
                volume_half1,
                volume_half2,
                apix=Apix,
                out_file=out_file,
                plot_file=plot_file,
            )
