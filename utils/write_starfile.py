'''
Create a Relion 3.0 star file from a particle stack and ctf parameters
'''

import argparse
import numpy as np
import sys, os
import pickle

import pandas as pd

from cryodrgn import dataset
from cryodrgn import utils
from cryodrgn import starfile
from cryodrgn import mrc
log = utils.log 

HEADERS = ['_rlnImageName',
           '_rlnDefocusU',
           '_rlnDefocusV',
           '_rlnDefocusAngle',
           '_rlnVoltage',
           '_rlnSphericalAberration',
           '_rlnAmplitudeContrast',
           '_rlnPhaseShift']

POSE_HDRS = ['_rlnAngleRot',
             '_rlnAngleTilt',
             '_rlnAnglePsi',
             '_rlnOriginX',
             '_rlnOriginY']

MICROGRAPH_HDRS = ['_rlnMicrographName',
                   '_rlnCoordinateX',
                   '_rlnCoordinateY']

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Input particles (.mrcs, .txt, .star, .cs)')
    parser.add_argument('ctf', help='Input ctf.pkl')
    parser.add_argument('--poses', help='Optionally include pose.pkl') 
    parser.add_argument('--ind', help='Optionally filter by selected index array (.pkl)')
    parser.add_argument('--datadir', type=os.path.abspath, help='Path prefix to particle stack if loading relative paths from a .star or .cs file')
    parser.add_argument('--full-path', action='store_true', help='Write the full path to particles (default: relative paths)')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output .star file')

    group = parser.add_argument_group('Optionally include additional star file columns')
    group.add_argument('--ref-star', help='Reference star file from original import')
    group.add_argument('--ref-star-relion31', action='store_true', help='Flag for relion3.1 star format')
    group.add_argument('--keep-micrograph', action='store_true', help='Include micrograph coordinate headers')
    return parser

def main(args):
    assert args.o.endswith('.star')
    particles = dataset.load_particles(args.particles, lazy=True, datadir=args.datadir)
    ctf = utils.load_pkl(args.ctf)
    assert ctf.shape[1] == 9, "Incorrect CTF pkl format"
    assert len(particles) == len(ctf), f"{len(particles)} != {len(ctf)}, Number of particles != number of CTF paraameters"
    if args.poses:
        poses = utils.load_pkl(args.poses)
        assert len(particles) == len(poses[0]), f"{len(particles)} != {len(poses)}, Number of particles != number of poses"
    log('{} particles'.format(len(particles)))

    if args.ref_star:
        ref_star = starfile.Starfile.load(args.ref_star, relion31=args.ref_star_relion31)
        assert len(ref_star) == len(particles), f"Particles in {args.particles} must match {args.ref_star}"

    if args.ind:
        ind = utils.load_pkl(args.ind)
        log(f'Filtering to {len(ind)} particles')
        particles = [particles[ii] for ii in ind]
        ctf = ctf[ind]
        if args.poses: 
            poses = (poses[0][ind], poses[1][ind])
        if args.ref_star:
            ref_star.df = ref_star.df.loc[ind] # note this does not reset the indexing 
    else:
        ind = np.arange(len(particles))

    ind += 1 # CHANGE TO 1-BASED INDEXING
    image_names = [img.fname for img in particles]
    if args.full_path:
        image_names = [os.path.abspath(img.fname) for img in particles]
    names = [f'{i}@{name}' for i,name in zip(ind, image_names)]

    ctf = ctf[:,2:]

    # convert poses
    if args.poses:
        eulers = utils.R_to_relion_scipy(poses[0]) 
        D = particles[0].get().shape[0]
        trans = poses[1] * D # convert from fraction to pixels

    data = {HEADERS[0]:names}
    for i in range(7):
        data[HEADERS[i+1]] = ctf[:,i]
    if args.poses:
        for i in range(3):
            data[POSE_HDRS[i]] = eulers[:,i]
        for i in range(2):
            data[POSE_HDRS[3+i]] = trans[:,i]
    df = pd.DataFrame(data=data) 
    headers = HEADERS + POSE_HDRS if args.poses else HEADERS
    if args.keep_micrograph:
        assert args.ref_star, "Must provide reference .star file with micrograph coordinates"
        log(f'Copying micrograph coordinates from {args.ref_star}')
        # TODO: Prepend path from args.ref_star to MicrographName?
        for h in MICROGRAPH_HDRS:
            df[h] = ref_star.df[h]
        headers += MICROGRAPH_HDRS

    s = starfile.Starfile(headers,df)
    s.write(args.o)

if __name__ == '__main__':
    main(parse_args().parse_args())
