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


def add_args(parser):
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
        help="Backproject the first N images (default: all images)",
    )
    parser.add_argument(
        "--log-interval",
        type=str,
        default="100",
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

@torch.inference_mode()
def _prepare_corner_indices_and_weights(ff_coord_T_3xN: torch.Tensor, D: int):
    """
    Vectorized version of the 8-corner 'splat' precomputation.

    ff_coord_T_3xN: (3, Npts) float coords in voxel-space, clipped to [-D/2, D/2].
    Returns:
      linear_idx_8N: (8, Npts) int64 flattened indices into a D^3 volume
      w_8N:          (8, Npts) float32 weights, w = 1 - ||dist||, clamped to [0, ∞)
    """
    d2 = int(D // 2)

    # floor/ceil for each axis, shape (N,)
    xf, yf, zf = ff_coord_T_3xN.floor().long()
    xc, yc, zc = ff_coord_T_3xN.ceil().long()

    # 8 corner combinations (xi, yi, zi), each shape (N,)
    xi = torch.stack([xf, xc, xf, xf, xc, xf, xc, xc], dim=0)
    yi = torch.stack([yf, yf, yc, yf, yc, yc, yf, yc], dim=0)
    zi = torch.stack([zf, zf, zf, zc, zf, zc, zc, zc], dim=0)

    # distances to each corner: sqrt((dx)^2+(dy)^2+(dz)^2), shapes (8, N)
    fx, fy, fz = ff_coord_T_3xN  # each (N,)
    dx = xi.to(torch.float32) - fx.unsqueeze(0)
    dy = yi.to(torch.float32) - fy.unsqueeze(0)
    dz = zi.to(torch.float32) - fz.unsqueeze(0)
    dist = (dx * dx + dy * dy + dz * dz).sqrt()

    w = (1.0 - dist).clamp_min_(0.0)  # match original w definition and clamp

    # convert (xi,yi,zi) to flattened linear indices for volume.view(-1)
    # target index = (zi + d2) * D*D + (yi + d2) * D + (xi + d2)
    xi = (xi + d2).clamp_(0, D - 1)
    yi = (yi + d2).clamp_(0, D - 1)
    zi = (zi + d2).clamp_(0, D - 1)

    linear_idx = zi * (D * D) + yi * D + xi
    return linear_idx.to(torch.int64), w.to(torch.float32)


@torch.inference_mode()
def _accumulate_slice(volume_1d: torch.Tensor,
                      counts_1d: torch.Tensor,
                      linear_idx_8N: torch.Tensor,
                      w_8N: torch.Tensor,
                      ff_eff_N: torch.Tensor,
                      ctfmul_sq_N: torch.Tensor):
    """
    Single scatter_add for volume and counts given precomputed corner indices/weights.

    volume_1d, counts_1d: D^3 flattened (float32) on device
    linear_idx_8N: (8, N), int64
    w_8N:          (8, N), float32
    ff_eff_N:           (N), float32  (ff * ctf_mul)
    ctfmul_sq_N:        (N), float32  (ctf_mul ** 2)
    """
    # Expand per-pixel scalars across the 8 corners
    vol_contrib_8N = w_8N * ff_eff_N.unsqueeze(0)     # (8, N)
    cnt_contrib_8N = w_8N * ctfmul_sq_N.unsqueeze(0)  # (8, N)

    # Flatten and scatter-add
    idx = linear_idx_8N.reshape(-1)
    volume_1d.scatter_add_(0, idx, vol_contrib_8N.reshape(-1))
    counts_1d.scatter_add_(0, idx, cnt_contrib_8N.reshape(-1))


def add_slice(volume, counts, ff_coord, ff, D, ctf_mul):
    d2 = int(D / 2)
    ff_coord = ff_coord.transpose(0, 1).clip(-d2, d2)
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


def regularize_volume(
    volume: torch.Tensor, counts: torch.Tensor, reg_weight: float
) -> torch.Tensor:
    regularized_counts = counts + reg_weight * counts.mean()
    regularized_counts *= counts.mean() / regularized_counts.mean()
    reg_volume = volume / regularized_counts

    return fft.ihtn_center(reg_volume[0:-1, 0:-1, 0:-1].cpu())


def main(args):
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

    # load the particles
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
            args.ntilts,
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
    Nimg = data.N
    lattice = Lattice(D, extent=D // 2, device=device)
    posetracker = PoseTracker.load(args.poses, Nimg, D, None, indices, device=device)

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
    img_iterator = range(min(args.first, Nimg)) if args.first else range(Nimg)

    # Keep masked coords on device once; shape (Npts, 3)
    coords_masked = lattice.coords[lattice_mask].to(device)

    if args.tilt:
        use_tilts = set(range(args.ntilts))
        img_iterator = [
            ii for ii in img_iterator if int(data.tilt_numbers[ii].item()) in use_tilts
        ]

    infer_ctx = torch.inference_mode()  # no-op on older torch if unavailable

    # Initialize tensors that will store backprojection results
    img_count = len(img_iterator)
    volume_full = torch.zeros((D, D, D), device=device, dtype=torch.float32)
    counts_full = torch.zeros((D, D, D), device=device, dtype=torch.float32)
    volume_half1 = torch.zeros((D, D, D), device=device, dtype=torch.float32)
    counts_half1 = torch.zeros((D, D, D), device=device, dtype=torch.float32)
    volume_half2 = torch.zeros((D, D, D), device=device, dtype=torch.float32)
    counts_half2 = torch.zeros((D, D, D), device=device, dtype=torch.float32)

    # Pre-flatten views to avoid repeated .view(-1) inside loop
    vol_full_1d = volume_full.view(-1)
    cnt_full_1d = counts_full.view(-1)
    if args.half_maps:
        vol_h1_1d = volume_half1.view(-1)
        cnt_h1_1d = counts_half1.view(-1)
        vol_h2_1d = volume_half2.view(-1)
        cnt_h2_1d = counts_half2.view(-1)


    # Figure out how often we are going to report progress w.r.t. images processed
    if args.log_interval == "auto":
        log_interval = max(round((img_count // 1000), -2), 100)
    elif args.log_interval.isnumeric():
        log_interval = int(args.log_interval)
    else:
        raise ValueError(
            f"Unrecognized argument for --log-interval: `{args.log_interval}`\n"
            f"Given value must be an integer or the label 'auto'!"
        )

    with infer_ctx:
        for i, ii in enumerate(img_iterator):
            if i % log_interval == 0:
                logger.info(f"fimage {ii} — {(i / img_count * 100):.1f}% done")

            # Pose and data fetch
            r, t = posetracker.get_pose(ii)  # r: (3,3), t: (2,) or None
            ff_tup = data.get_tilt(ii) if args.tilt else data[ii]
            assert isinstance(ff_tup, tuple)
            ff = ff_tup[0].to(device)  # (D,D)
            ff = ff.view(-1)[lattice_mask]  # (Npts,)

            # CTF / filters (per-pixel)
            ctf_mul = 1.0
            if ctf_params is not None:
                freqs = lattice.freqs2d / ctf_params[ii, 0]
                c = ctf.compute_ctf(freqs, *ctf_params[ii, 1:]).view(-1)[lattice_mask]
                if args.ctf_alg == "flip":
                    ff = ff * c.sign()
                else:
                    ctf_mul = c

            # In-plane real-space translation (phase shift in Fourier domain)
            if t is not None:
                ff = lattice.translate_ht(ff.view(1, -1), t.view(1, 1, 2), lattice_mask).view(-1)

            # Dose filters for tilt data
            if args.tilt:
                # build once per ii without creating large temporaries
                tilt_idxs = torch.tensor([ii], device=device)
                dose_filters = data.get_dose_filters(tilt_idxs, lattice, ctf_params[ii, 0])[0]
                ctf_mul = ctf_mul * dose_filters[lattice_mask]

            # Precompute 8-corner indices & weights once for this image
            # coords_masked: (Npts, 3). We need (3, Npts) for the helper.
            ff_coord_T = (coords_masked @ r).T.clamp_(min=-D/2, max=D/2)  # (3, Npts)
            linear_idx_8N, w_8N = _prepare_corner_indices_and_weights(ff_coord_T, D)

            # Expand per-pixel scalars; reuse for full + half
            ff_eff    = (ff * ctf_mul).to(torch.float32)      # (Npts,)
            ctfmul_sq = (ctf_mul * ctf_mul).to(torch.float32) # (Npts,)

            # Accumulate into full map
            _accumulate_slice(vol_full_1d, cnt_full_1d, linear_idx_8N, w_8N, ff_eff, ctfmul_sq)

            # Accumulate into the appropriate half-map without recomputing indices/weights
            if args.half_maps:
                if (ii % 2) == 0:
                    _accumulate_slice(vol_h1_1d, cnt_h1_1d, linear_idx_8N, w_8N, ff_eff, ctfmul_sq)
                else:
                    _accumulate_slice(vol_h2_1d, cnt_h2_1d, linear_idx_8N, w_8N, ff_eff, ctfmul_sq)

    td = time.time() - t1
    logger.info(
        f"Backprojected {img_count} images "
        f"in {td:.2f}s ({(td / img_count):4f}s per image)"
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
