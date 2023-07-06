"""Classify a set of maps by comparing to known volumes"""

import argparse
import logging
import numpy as np
from cryodrgn import fft
from cryodrgn.source import ImageSource
import torch

logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vols1", help="Input")
    parser.add_argument("vols2", help="Input")
    parser.add_argument("--mask")
    parser.add_argument("--do-mse", action="store_true")
    # parser.add_argument("-o", help="Output")
    return parser


# TODO: implement diffmap normalization method here
# Currently just rescales values to between 10th and 90th percentiles
def normalize_vols(vols, vol_ref):
    quartiles = torch.quartile(vol_ref, torch.tensor([0.1, 0.9]))
    new_range = quartiles[1] - quartiles[0]
    new_min = quartiles[0]
    for ii, vol in vols:
        min_val = torch.min(vol)
        max_val = torch.max(vol)
        vols[ii,:] = new_range * (vol - min_val)/(max_val - min_val) + new_min

    return vols


def get_fsc(vol1, vol2, r, cutoff=cutoff):
    vol1 = fft.fftn_center(vol1)
    vol2 = fft.fftn_center(vol2)

    prev_mask = np.zeros((D, D, D), dtype=bool)
    fsc = [1.0]
    for i in range(1, D // 2):
        mask = r < i
        shell = np.where(mask & np.logical_not(prev_mask))
        v1 = vol1[shell]
        v2 = vol2[shell]
        p = np.vdot(v1, v2) / (np.vdot(v1, v1) * np.vdot(v2, v2)) ** 0.5
        fsc.append(float(p.real))
        prev_mask = mask
    fsc = np.asarray(fsc)
    x = np.arange(D // 2) / D

    res = 1.0/x[np.where(fsc < cutoff)[0][0]]
    return res

def get_pairwise_fsc(vols1, vols2, cutoff=0.5):
    x = np.arange(-D // 2, D // 2)
    x2, x1, x0 = np.meshgrid(x, x, x, indexing="ij")
    coords = np.stack((x0, x1, x2), -1)
    r = (coords**2).sum(-1) ** 0.5

    fsc = np.array(vols1.shape[0], vols2.shape[0])
    for ii, vol1 in enumerate(vols1):
        for jj, vol2 in enumerate(vols2):
            fsc[ii, jj] = get_fsc(vol1, vol2, r, cutoff=cutoff)

    return torch.from_numpy(fsc) 

def main(args):
    vols1 = ImageSource.from_file(args.vols1)
    vols2 = ImageSource.from_file(args.vols2)
    
    D = vols1.D
    assert D == vols2.D

    vols1 = vols1.images().view(-1, D*D*D)
    vols2 = vols2.images().view(-1, D*D*D)


    if args.mask:
        mask = ImageSource.from_file(args.mask)
        assert D == mask.D

        mask = mask.images()
        mask = mask.to(torch.bool)

        vols1 *= mask.repeat(vols1.shape[0], 1, 1, 1)
        vols2 *= mask.repeat(vols2.shape[0], 1, 1, 1).view(-1, D*D*D)

    vol_ref = vols2[0:]
    vols1_norm = normalize_vols(vols1, vol_ref).view(-1, D*D*D)
    vols2_norm = normalize_vols(vols2, vol_ref).view(-1, D*D*D)
    mse = torch.cdist(vols1_norm, vols2_norm, p=2)**2/(D*D*D)
    
    fsc = get_pairwise_fsc(vols1, vols2, D)

    classes = torch.argmax(fsc, dim=1)
    if args.do_mse:
        classes = torch.argmax(mse, dim=1)

    num_classes = []
    for ii in range(vols2.shape[0]):
        num_classes += [sum(classes == ii)]

    print(num_classes)

if __name__ == "__main__":
    main(parse_args().parse_args())

