"""
Create a Relion 3.0 star file from a particle stack and ctf parameters
"""

import argparse
import os, sys
import numpy as np
import pandas as pd
import random
import logging
from cryodrgn import utils
from cryodrgn.source import ImageSource, StarfileSource
from cryodrgn.starfile import Starfile
import pickle
import re
# file="/scratch/gpfs/ZHONGE/mj7341/research/00_moml/antibody/dataset/conformational"
# python /scratch/gpfs/ZHONGE/mj7341/MoML/cryosim/cryodrgn/cryodrgn/commands_utils/write_star.py $file/add_noise/128_chimera_resample/snr01/mrcs/sorted_particles.128.txt --ctf $file/integrated_ctf.pkl --poses $file/integrated_poses_chimera.pkl -o $file/add_noise/128_chimera_resample/snr01/snr01_star.star 
log = print
logger = logging.getLogger(__name__)


CTF_HEADERS = [
    "_rlnDefocusU",
    "_rlnDefocusV",
    "_rlnDefocusAngle",
    "_rlnVoltage",
    "_rlnSphericalAberration",
    "_rlnAmplitudeContrast",
    "_rlnPhaseShift",
]

POSE_HDRS = [
    "_rlnAngleRot",
    "_rlnAngleTilt",
    "_rlnAnglePsi",
    "_rlnOriginX",
    "_rlnOriginY",
]

RELION_HEADERS = ["_rlnOpticsGroupName","_rlnOpticsGroup", "_rlnMicrographOriginalPixelSize", "_rlnVoltage",
            "_rlnSphericalAberration", "_rlnAmplitudeContrast", "_rlnImagePixelSize", "_rlnImageSize", "_rlnImageDimensionality",
            "_rlnCtfDataAreCtfPremultiplied"]

CRYOSPARC_HEADERS = ['_rlnImagePixelSize', '_rlnImageSize', '_rlnImageDimensionality', '_rlnOpticsGroup', '_rlnRandomSubset']

def add_args(parser):
    parser.add_argument("particles", help="Input particles (.mrcs, .txt, .star)")
    parser.add_argument("--ctf", help="Input ctf.pkl")
    parser.add_argument("--poses", help="Optionally include pose.pkl")
    parser.add_argument(
        "--ind", help="Optionally filter by array of selected indices (.pkl)"
    )
    parser.add_argument(
        "--full-path",
        action="store_true",
        help="Write the full path to particles (default: relative paths)",
    )
    parser.add_argument(
        "-o", type=os.path.abspath, required=True, help="Output .star file"
    )
    parser.add_argument(
        "--img_dim",
        type=int,
        default=2,
        help="_rlnImageDimensionality (default: %(default)s)",
    )
    parser.add_argument(
        "--optics_group",
        type=int,
        default=1,
        help="_rlnOpticsGroup (default: %(default)s)",
    )
    
    return parser

def natural_sort_key(s):
    # Convert the string to a list of text and numbers
    parts = re.split('([0-9]+)', s)
    
    # Convert numeric parts to integers for proper numeric comparison
    parts[1::2] = map(int, parts[1::2])
    
    return parts

# def print_ctf_params(params: np.ndarray) -> None:
#     assert len(params) == 9
#     logger.info("Image size (pix)  : {}".format(int(params[0])))
#     logger.info("A/pix             : {}".format(params[1]))
#     logger.info("DefocusU (A)      : {}".format(params[2]))
#     logger.info("DefocusV (A)      : {}".format(params[3]))
#     logger.info("Dfang (deg)       : {}".format(params[4]))
#     logger.info("voltage (kV)      : {}".format(params[5]))
#     logger.info("cs (mm)           : {}".format(params[6]))
#     logger.info("w                 : {}".format(params[7]))
#     logger.info("Phase shift (deg) : {}".format(params[8]))

def main(args):
    assert args.o.endswith(".star"), "Output file must be .star file"
    input_ext = os.path.splitext(args.particles)[-1]
    assert input_ext in (
        ".mrcs",
        ".txt",
        ".star",
    ), "Input file must be .mrcs/.txt/.star"

    # Either accept an input star file, or an input .mrcs/.txt with optional ctf/pose pkl file(s)
    ctf = poses = eulers = trans = None
    if input_ext == ".star":
        assert (
            args.poses is None
        ), "--poses cannot be specified when input is a starfile (poses are obtained from starfile)"
        assert (
            args.ctf is None
        ), "--ctf cannot be specified when input is a starfile (ctf information are obtained from starfile)"

    particles = ImageSource.from_file(args.particles, lazy=True)

    if args.ctf:
        ctf = utils.load_pkl(args.ctf)
        assert ctf.shape[1] == 9, "Incorrect CTF pkl format"
        assert len(particles) == len(
            ctf
        ), f"{len(particles)} != {len(ctf)}, Number of particles != number of CTF parameters"
    if args.poses:
        poses = utils.load_pkl(args.poses)
        assert len(particles) == len(
            poses[0]
        ), f"{len(particles)} != {len(poses)}, Number of particles != number of poses"
    logger.info(f"{len(particles)} particles in {args.particles}")

    ind = np.arange(particles.n)
    if args.ind:
        ind = utils.load_pkl(args.ind)
        logger.info(f"Filtering to {len(ind)} particles")
        if ctf is not None:
            ctf = ctf[ind]
        if poses is not None:
            poses = (poses[0][ind], poses[1][ind])

    if input_ext == ".star":
        assert isinstance(particles, StarfileSource)
        df = particles.df.loc[ind]
    else:
        image_names = particles.filenames[ind]
        particle_path, _ = os.path.split(args.particles)
        file_names_lst = os.listdir(particle_path)
        file_names_lst = sorted(file_names_lst, key=natural_sort_key)
        if args.full_path:
            image_names = [os.path.abspath(image_name) for image_name in image_names]
        names = []
        j=1
        for i, name in zip(ind, image_names):
            if j % 1000 ==1:
                j=1
            print(f"{j}@{name}"+particle_path+'/'+file_names_lst[i//1000])
            if file_names_lst[i//1000].split('.')[-1] == 'txt':
                continue
            else:
                names.append(f"{j}@{name}"+particle_path+'/'+file_names_lst[i//1000])
            j = j+1

        # convert poses
        if poses is not None:
            eulers = utils.R_to_relion_scipy(poses[0])
            D = particles[0].shape[-1]
            trans = poses[1] * D  # convert from fraction to pixels
        
        # Create a new dataframe with required star file headers
        data = {"_rlnImageName": names}
        if ctf is not None:
            for i in range(7):
                data[CTF_HEADERS[i]] = ctf[:, i+2]

        if eulers is not None and trans is not None:
            for i in range(3):
                data[POSE_HDRS[i]] = eulers[:, i]  # type: ignore
            for i in range(2):
                data[POSE_HDRS[3 + i]] = trans[:, i]
        
        #### CryoSPARC ####
        column_values = [ctf[0][1], ctf[0][0], args.img_dim, args.optics_group, 1]
        for i in range(len(CRYOSPARC_HEADERS)):
            data[CRYOSPARC_HEADERS[i]] = column_values[i]
        df = pd.DataFrame(data=data)
        # Half Set
        lst = list(range(len(df)))
        list_sampled = random.sample(lst, k=len(lst)//2)
        df[CRYOSPARC_HEADERS[-1]].loc[list_sampled] = 2
    s = Starfile(headers=None, df=df)

    #### Relion ####
    relion_values = ["opticsGroup1", args.optics_group, ctf[0][1], ctf[0][5], ctf[0][6],
                    ctf[0][7], ctf[0][1], ctf[0][0], args.img_dim, 0]
    s.write(args.o, RELION_HEADERS, relion_values)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
