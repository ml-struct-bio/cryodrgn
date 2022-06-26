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

def add_args(parser):
    parser.add_argument('particles', help='Input particles (.mrcs, .txt)')
    parser.add_argument('ctf', help='Input ctf.pkl')
    parser.add_argument('--poses', help='Optionally include pose.pkl') 
    parser.add_argument('--ind', help='Optionally filter by array of selected indices (.pkl)')
    parser.add_argument('--full-path', action='store_true', help='Write the full path to particles (default: relative paths)')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output .star file')

    group = parser.add_argument_group('Optionally include additional star file columns')
    group.add_argument('--ref-star', help='Reference star file from original import')
    group.add_argument('--keep-micrograph', action='store_true', help='Include micrograph coordinate headers')
    group.add_argument('--copy-header', nargs='*', help='Additional headers to copy from the reference star file')
    return parser

def parse_chunk_size(txtfile):
    lines = open(txtfile,'r').readlines()
    def abspath(f):
        if os.path.isabs(f):
            return f
        base = os.path.dirname(os.path.abspath(txtfile))
        return os.path.join(base,f)
    lines = [abspath(x).strip() for x in lines]
    return [mrc.parse_header(f).fields['nz'] for f in lines]

def main(args):
    assert args.o.endswith('.star'), "Output file must be .star file"
    assert args.particles.endswith('.mrcs') or args.particles.endswith('.txt'), "Input file must be .mrcs or .txt"

    particles = dataset.load_particles(args.particles, lazy=True)
    ctf = utils.load_pkl(args.ctf)
    assert ctf.shape[1] == 9, "Incorrect CTF pkl format"
    assert len(particles) == len(ctf), f"{len(particles)} != {len(ctf)}, Number of particles != number of CTF paraameters"
    if args.poses:
        poses = utils.load_pkl(args.poses)
        assert len(particles) == len(poses[0]), f"{len(particles)} != {len(poses)}, Number of particles != number of poses"
    log(f'{len(particles)} particles in {args.particles}')

    if args.ref_star:
        ref_star = starfile.Starfile.load(args.ref_star)
        assert len(ref_star) == len(particles), f"{len(particles)} != {len(ref_star)}, Number of particles in {args.particles} != number of particles in {args.ref_star}"

    # Get index for particles in each .mrcs file
    if args.particles.endswith('.txt'):
        N_per_chunk = parse_chunk_size(args.particles)
        particle_ind = np.concatenate([np.arange(nn) for nn in N_per_chunk])
        assert len(particle_ind) == len(particles)
    else: # single .mrcs file
        particle_ind = np.arange(len(particles))

    if args.ind:
        ind = utils.load_pkl(args.ind)
        log(f'Filtering to {len(ind)} particles')
        particles = [particles[ii] for ii in ind]
        ctf = ctf[ind]
        if args.poses: 
            poses = (poses[0][ind], poses[1][ind])
        if args.ref_star:
            ref_star.df = ref_star.df.loc[ind] 
            # reset the index in the dataframe to avoid any downstream indexing issues
            ref_star.df.reset_index(inplace=True)
        particle_ind = particle_ind[ind]

    particle_ind += 1 # CHANGE TO 1-BASED INDEXING
    image_names = [img.fname for img in particles]
    if args.full_path:
        image_names = [os.path.abspath(img.fname) for img in particles]
    names = [f'{i}@{name}' for i,name in zip(particle_ind, image_names)]

    ctf = ctf[:,2:]

    # convert poses
    if args.poses:
        eulers = utils.R_to_relion_scipy(poses[0]) 
        D = particles[0].get().shape[0]
        trans = poses[1] * D # convert from fraction to pixels

    # Create a new dataframe with required star file headers
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
            log(f'  Copying {h}')
            df[h] = ref_star.df[h]
        headers += MICROGRAPH_HDRS

    if args.copy_header is not None:
        assert args.ref_star, "Must provide reference .star file"
        log(f'Copying additional columns from {args.ref_star}')
        for h in args.copy_header:
            log(f'  Copying {h}')
            df[h] = ref_star.df[h]
            headers.append(h)

    s = starfile.Starfile(headers,df)
    s.write(args.o)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
