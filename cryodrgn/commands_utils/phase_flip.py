"""Phase flip images by CTF sign"""

import argparse
import os
import logging
import numpy as np
import torch
from cryodrgn import ctf, fft
from cryodrgn.mrc import MRCFile
from cryodrgn.source import ImageSource
from cryodrgn.utils import meshgrid_2d

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument("mrcs", help="Input particles (.mrcs, .txt, .star, or .cs)")
    parser.add_argument("ctf_params", help="Input CTF parameters (.pkl)")
    parser.add_argument(
        "--datadir",
        help="Optionally overwrite path to starfile .mrcs if loading from a starfile",
    )
    parser.add_argument("-o", type=os.path.abspath, help="Output .mrcs")
    return parser


def main(args):
    imgs = ImageSource.from_file(args.mrcs, lazy=True, datadir=args.datadir)
    D = imgs.D
    ctf_params = ctf.load_ctf_for_training(D, args.ctf_params)
    ctf_params = torch.Tensor(ctf_params)
    assert len(imgs) == len(ctf_params), f"{len(imgs)} != {len(ctf_params)}"

    fx2, fy2 = meshgrid_2d(-0.5, 0.5, D, endpoint=False)
    freqs = torch.stack((fx2.ravel(), fy2.ravel()), dim=1)

    def compute_ctf(freqs, ctf_params):
        dfu = ctf_params[:, [0]]
        dfv = ctf_params[:, [1]]
        dfang = ctf_params[:, [2]]
        volt = ctf_params[:, [3]]
        cs = ctf_params[:, [4]]
        w = ctf_params[:, [5]]
        phase_shift = ctf_params[:, [6]]

        bfactor = None
        if ctf_params.shape[1] >= 8:
            bfactor = ctf_params[:, [7]]

        return ctf.compute_ctf(
            freqs, dfu, dfv, dfang, volt, cs, w, phase_shift, bfactor
        )

    def transform_fn(img, indices):
        _freqs = freqs.expand(len(indices), -1, -1)
        _ctf_params = ctf_params[..., np.newaxis, np.newaxis]
        c = compute_ctf(_freqs / _ctf_params[indices, 0], ctf_params[indices, 1:])
        c = c.reshape((-1, D, D))
        ff = fft.fft2_center(img)
        ff *= torch.sign(c)
        img2 = fft.ifftn_center(ff)

        return img2.cpu().numpy().astype(np.float32)

    logger.info(f"Writing {args.o}")
    MRCFile.write(args.o, imgs, transform_fn=transform_fn)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
