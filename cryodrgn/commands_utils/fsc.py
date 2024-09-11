"""Compute Fourier shell correlations between two volumes, applying an optional mask.

When using `--ref-volume`, this reference volume will be used to create a loose and a
tight mask using dilation and cosine edges; the latter mask will be corrected using
phase randomization as implemented in cryoSPARC.

See also
--------
`cryodrgn.commands_utils.plot_fsc` — for just plotting already-calculated FSCs
`cryodrgn.commands_utils.gen_mask` — for generating custom dilation + cosine edge masks

Example usage
-------------
$ cryodrgn_utils fsc volume1.mrc volume2.mrc

# Save FSC values to file and produce an FSC plot
$ cryodrgn_utils fsc vol1.mrc vol2.mrc -o fsc.txt -p

# Also apply a mask before computing FSCs, and also produce a plot of FSC curves.
$ cryodrgn_utils fsc vol1.mrc vol2.mrc --mask test-mask.mrc -o fsc.txt -p

# Also apply phase randomization at Fourier shells for resolutions < 10 angstroms
$ cryodrgn_utils fsc vol1.mrc vol2.mrc --mask test-mask.mrc -o fsc.txt -p fsc-plot.png
                                       --corrected 10

# Do cryoSPARC-style phase randomization with a tight mask; create an FSC plot with
# curves for no mask, spherical mask, loose mask, tight mask, and corrected tight mask,
# with loose and tight masks generated using the (first-given) full volume.
$ cryodrgn_utils fsc vol1.mrc vol2.mrc --ref-volume=full.mrc -p fsc-plot.png

# Do phase randomization as above using the given mask instead of the default tight
# mask, which will also be used to calculate corrected FSCs.
$ cryodrgn_utils fsc half_vol_a.mrc half_vol_b.mrc --ref-volume=backproject.mrc \
                     -p tighter-mask.png mask=tighter-mask.mrc

"""
import os
import argparse
import logging
import numpy as np
import pandas as pd
import torch
from typing import Optional, Union
from cryodrgn import fft
from cryodrgn.source import MRCFileSource
from cryodrgn.commands_utils.plot_fsc import create_fsc_plot
from cryodrgn.masking import spherical_window_mask, cosine_dilation_mask

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments available for use with `cryodrgn_utils fsc`."""
    parser.add_argument("volumes", nargs=2, help="the two .mrc volume files to compare")

    parser.add_argument(
        "--ref-volume",
        type=os.path.abspath,
        help="if given, create a cryoSPARC-style FSC plot instead using this .mrc "
        "volume as the reference to produce masks",
    )
    parser.add_argument(
        "--mask",
        type=os.path.abspath,
        help="if given, apply the mask in this file before calculating half-map FSCs, "
        "or use it for the correction step if calculating phase-randomized FSCs",
    )
    parser.add_argument(
        "--corrected",
        type=float,
        help="use cryoSPARC-style high resolution phase randomization beyond this "
        "resolution (in angstroms) to correct for possible effects of tight masking",
    )
    parser.add_argument(
        "--plot",
        "-p",
        type=os.path.abspath,
        nargs="?",
        const=True,
        default=None,
        help="also plot the FSC curve; optionally supply a .png file name instead of "
        "generating one automatically",
    )
    parser.add_argument(
        "--Apix",
        type=float,
        help="Ang/pixels to use when printing the resolutions at thresholds to "
        "override those found in the volumes, or replace them if missing in volumes",
    )
    parser.add_argument(
        "--outtxt",
        "-o",
        type=os.path.abspath,
        help=(
            "if given, a file to save the FSC values, "
            "with each space-delimited row as <resolution> <fsc_val>; "
            "otherwise print these to screen"
        ),
    )


def get_fftn_center_dists(box_size: int) -> np.array:
    """Get distances from the center (and hence the resolution) for FFT co-ordinates."""

    x = np.arange(-box_size // 2, box_size // 2)
    x2, x1, x0 = np.meshgrid(x, x, x, indexing="ij")
    coords = np.stack((x0, x1, x2), -1)
    dists = (coords**2).sum(-1) ** 0.5
    assert dists[box_size // 2, box_size // 2, box_size // 2] == 0.0

    return dists


def calculate_fsc(
    v1: Union[np.ndarray, torch.Tensor], v2: Union[np.ndarray, torch.Tensor]
) -> float:
    """Calculate the Fourier Shell Correlation between two complex vectors."""
    var = (np.vdot(v1, v1) * np.vdot(v2, v2)) ** 0.5

    if var:
        fsc = float((np.vdot(v1, v2) / var).real)
    else:
        fsc = 1.0

    return fsc


def get_fsc_curve(
    vol1: torch.Tensor,
    vol2: torch.Tensor,
    initial_mask: Optional[torch.Tensor] = None,
    out_file: Optional[str] = None,
) -> pd.DataFrame:
    """Calculate the FSCs between two volumes across all available resolutions."""

    # Apply the given mask before applying the Fourier transform
    maskvol1 = vol1 * initial_mask if initial_mask is not None else vol1.clone()
    maskvol2 = vol2 * initial_mask if initial_mask is not None else vol2.clone()
    box_size = vol1.shape[0]
    dists = get_fftn_center_dists(box_size)
    maskvol1 = fft.fftn_center(maskvol1)
    maskvol2 = fft.fftn_center(maskvol2)

    prev_mask = np.zeros((box_size, box_size, box_size), dtype=bool)
    fsc = [1.0]
    for i in range(1, box_size // 2):
        mask = dists < i
        shell = np.where(mask & np.logical_not(prev_mask))
        fsc.append(calculate_fsc(maskvol1[shell], maskvol2[shell]))
        prev_mask = mask

    fsc_vals = pd.DataFrame(
        dict(pixres=np.arange(box_size // 2) / box_size, fsc=fsc), dtype=float
    )
    if out_file is not None:
        logger.info(f"Saving FSC values to {out_file}")
        fsc_vals.to_csv(out_file, sep=" ", header=True, index=False)

    return fsc_vals


def get_fsc_thresholds(
    fsc_vals: pd.DataFrame, apix: float, verbose: bool = True
) -> tuple[float, float]:
    """Retrieve the max resolutions at which an FSC curve is above 0.5 and 0.143."""

    if ((fsc_vals.pixres > 0) & (fsc_vals.fsc >= 0.5)).any():
        res_05 = fsc_vals.pixres[fsc_vals.fsc >= 0.5].max()
        if verbose:
            logger.info("res @ FSC=0.5: {:.4g} ang".format((1 / res_05) * apix))
    else:
        res_05 = None
        if verbose:
            logger.warning("res @ FSC=0.5: N/A")

    if ((fsc_vals.pixres > 0) & (fsc_vals.fsc >= 0.143)).any():
        res_143 = fsc_vals.pixres[fsc_vals.fsc >= 0.143].max()
        if verbose:
            logger.info("res @ FSC=0.143: {:.4g} ang".format((1 / res_143) * apix))
    else:
        res_143 = None
        if verbose:
            logger.warning("res @ FSC=0.143: N/A")

    return res_05, res_143


def randomize_phase(cval: complex) -> complex:
    """Create a new complex value with the same amplitude but scrambled phase."""
    amp = (cval.real**2.0 + cval.imag**2.0) ** 0.5
    angrand = np.random.random() * 2 * np.pi

    return complex(amp * np.cos(angrand), amp * np.sin(angrand))


def correct_fsc(
    fsc_vals: pd.DataFrame,
    vol1: torch.Tensor,
    vol2: torch.Tensor,
    randomization_threshold: float,
    initial_mask: Optional[torch.Tensor] = None,
) -> pd.DataFrame:
    """Apply phase-randomization null correction to given FSC volumes past a resolution.

    This function implements cryoSPARC-style correction to an FSC curve to account for
    the boost in FSCs that can be attributed to a mask that is too sharp or too tightly
    fits the volumes and thus introduces an artificial source of correlation.

    """
    box_size = vol1.shape[0]
    if fsc_vals.shape[0] != (box_size // 2):
        raise ValueError(
            f"Given FSC values must have (D // 2) + 1 = {(box_size // 2) + 1} entries, "
            f"instead have {fsc_vals.shape[0]}!"
        )

    # Randomize phases in the raw half-maps beyond the given threshold
    dists = get_fftn_center_dists(box_size)
    fftvol1 = fft.fftn_center(vol1)
    fftvol2 = fft.fftn_center(vol2)
    phase_res = int(randomization_threshold * box_size)
    rand_shell = np.where(dists >= phase_res)
    fftvol1[rand_shell] = fftvol1[rand_shell].apply_(randomize_phase)
    fftvol2[rand_shell] = fftvol2[rand_shell].apply_(randomize_phase)
    fftvol1 = fft.ifftn_center(fftvol1)
    fftvol2 = fft.ifftn_center(fftvol2)

    # Apply the given masks then go back into Fourier space
    maskvol1 = fftvol1 * initial_mask if initial_mask is not None else fftvol1.clone()
    maskvol2 = fftvol2 * initial_mask if initial_mask is not None else fftvol2.clone()
    maskvol1 = fft.fftn_center(maskvol1)
    maskvol2 = fft.fftn_center(maskvol2)

    # re-calculate the FSCs past the resolution using the phase-randomized volumes
    prev_mask = np.zeros((box_size, box_size, box_size), dtype=bool)
    fsc = fsc_vals.fsc.tolist()
    for i in range(1, box_size // 2):
        mask = dists < i
        shell = np.where(mask & np.logical_not(prev_mask))

        if i > phase_res:
            p = calculate_fsc(maskvol1[shell], maskvol2[shell])

            # normalize the original FSC value using the phase-randomized value
            if p == 1.0:
                fsc[i] = 0.0
            elif not np.isnan(p):
                fsc[i] = np.clip((fsc[i] - p) / (1 - p), 0, 1.0)

        prev_mask = mask

    return pd.DataFrame(
        dict(pixres=np.arange(box_size // 2) / box_size, fsc=fsc), dtype=float
    )


def calculate_cryosparc_fscs(
    full_vol: torch.Tensor,
    half_vol1: torch.Tensor,
    half_vol2: torch.Tensor,
    sphere_mask: Optional[Union[np.ndarray, torch.Tensor]] = None,
    loose_mask: tuple[int, int] = (25, 15),
    tight_mask: Union[tuple[int, int], np.ndarray] = (6, 6),
    apix: float = 1.0,
    out_file: Optional[str] = None,
    plot_file: Optional[str] = None,
) -> pd.DataFrame:
    """Calculating cryoSPARC-style FSC curves with phase randomization correction."""
    if sphere_mask is None:
        sphere_mask = spherical_window_mask(D=full_vol.shape[0])

    masks = {
        "No Mask": None,
        "Spherical": sphere_mask,
        "Loose": cosine_dilation_mask(
            full_vol, dilation=loose_mask[0], edge_dist=loose_mask[1], apix=apix
        ),
    }
    if isinstance(tight_mask, tuple):
        masks["Tight"] = cosine_dilation_mask(
            full_vol, dilation=tight_mask[0], edge_dist=tight_mask[1], apix=apix
        )
    elif isinstance(tight_mask, (np.ndarray, torch.Tensor)):
        masks["Tight"] = tight_mask
    else:
        raise TypeError(
            f"`tight_mask` must be an array or a tuple giving dilation and cosine edge "
            f"size in pixels, instead given {type(tight_mask).__name__}!"
        )

    fsc_vals = {
        mask_lbl: get_fsc_curve(half_vol1, half_vol2, initial_mask=mask)
        for mask_lbl, mask in masks.items()
    }
    fsc_thresh = {
        mask_lbl: get_fsc_thresholds(fsc_df, apix, verbose=False)[1]
        for mask_lbl, fsc_df in fsc_vals.items()
    }

    if fsc_thresh["Tight"] is not None:
        fsc_vals["Corrected"] = correct_fsc(
            fsc_vals["Tight"],
            half_vol1,
            half_vol2,
            randomization_threshold=0.75 * fsc_thresh["Tight"],
            initial_mask=masks["Tight"],
        )
        fsc_thresh["Corrected"] = get_fsc_thresholds(
            fsc_vals["Corrected"], apix, verbose=False
        )[1]
    else:
        fsc_vals["Corrected"] = fsc_vals["Tight"]
        fsc_thresh["Corrected"] = fsc_thresh["Tight"]

    # Report corrected FSCs by printing FSC=0.5 and FSC=0.143 threshold values to screen
    get_fsc_thresholds(fsc_vals["Corrected"], apix)

    if plot_file is not None:
        fsc_angs = {
            mask_lbl: ((1 / fsc_val) * apix) for mask_lbl, fsc_val in fsc_thresh.items()
        }
        fsc_plot_vals = {
            f"{mask_lbl}  ({fsc_angs[mask_lbl]:.2f}Å)": fsc_df
            for mask_lbl, fsc_df in fsc_vals.items()
        }
        create_fsc_plot(fsc_vals=fsc_plot_vals, outfile=plot_file, apix=apix)

    pixres_index = {tuple(vals.pixres.values) for vals in fsc_vals.values()}
    assert len(pixres_index) == 1
    pixres_index = tuple(pixres_index)[0]

    fsc_vals = pd.DataFrame(
        {k: vals.fsc.values for k, vals in fsc_vals.items()}, index=list(pixres_index)
    )
    fsc_vals.index.name = "pixres"

    if out_file is not None:
        fsc_vals.reset_index(inplace=True, drop=False)
        logger.info(f"Saving FSC values to {out_file}")
        fsc_vals.columns = fsc_vals.columns.str.replace(" ", "")
        fsc_vals.round(6).to_csv(out_file, sep=" ", header=True, index=False)

    return fsc_vals


def main(args: argparse.Namespace) -> None:
    """Calculate FSC curves based on command-line arguments (see `add_args()` above)."""

    volumes = [MRCFileSource(vol_file) for vol_file in args.volumes]
    mask = MRCFileSource(args.mask).images() if args.mask is not None else None

    if args.Apix:
        apix = args.Apix
    else:
        vol_apixs = {vol.apix for vol in volumes if vol.apix is not None}

        if len(vol_apixs) == 0:
            apix = 1.0
        elif len(vol_apixs) == 1:
            apix = tuple(vol_apixs)[0]
        else:
            raise ValueError(
                f"These volumes have different A/px values: {vol_apixs}"
                f"\nUse `--Apix` to supply your own A/px value to override this!"
            )

    if args.plot:
        if isinstance(args.plot, bool):
            if args.outtxt:
                plot_file = "".join([os.path.splitext(args.outtxt)[0], ".png"])
            else:
                plot_file = "fsc-plot.png"
        else:
            plot_file = str(args.plot)
    else:
        plot_file = None

    logger.info(
        f"Calculating FSCs between `{args.volumes[0]}` and `{args.volumes[1]}`..."
    )
    if args.ref_volume is None:
        fsc_vals = get_fsc_curve(
            volumes[0].images(),
            volumes[1].images(),
            mask,
            out_file=args.outtxt,
        )
        _ = get_fsc_thresholds(fsc_vals, apix)

        if args.corrected is not None:
            if args.corrected >= 1:
                args.corrected = (args.corrected / apix) ** -1

            fsc_vals = correct_fsc(
                fsc_vals, volumes[0].images(), volumes[1].images(), args.corrected, mask
            )
        if args.plot:
            create_fsc_plot(fsc_vals=fsc_vals, outfile=plot_file, apix=apix)

    else:
        logger.info(
            f"Using `{args.ref_volume}` as the reference volume for creating masks..."
        )
        if args.corrected is not None:
            raise ValueError(
                "Cannot provide your own phase randomization threshold as this will "
                "be computed for you when using a `--ref-volume`!"
            )

        ref_volume = MRCFileSource(args.ref_volume)
        if mask is None:
            mask = (6, 6)

        fsc_vals = calculate_cryosparc_fscs(
            ref_volume.images(),
            volumes[0].images(),
            volumes[1].images(),
            tight_mask=mask,
            apix=apix,
            out_file=args.outtxt,
            plot_file=plot_file,
        )

    # If we didn't save the FSC values to file, print them to screen instead
    if not args.outtxt:
        show_indx = "pixres" not in fsc_vals.columns
        if show_indx:
            fsc_vals.index = np.round(fsc_vals.index, 3)

        fsc_str = fsc_vals.round(4).to_csv(sep="\t", index=show_indx, header=True)
        logger.info(f"\n{fsc_str}")
