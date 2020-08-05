'''
Create a Relion 3.0 star file from a particle stack and ctf parameters
'''

import argparse
import numpy as np
import sys, os
import pickle

import pandas as pd

from cryodrgn import utils
from cryodrgn import starfile
from cryodrgn import mrc
log = utils.log 

HEADERS = ['_rlnImageName',
           '_rlnDefocusU',
           '_rlnDefocusV',
           '_rlnDefocusAngle',
           '_rlnVoltage',
           '_rlnAmplitudeContrast',
           '_rlnSphericalAberration',
           '_rlnPhaseShift']

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', type=os.path.abspath, help='Input .mrcs')
    parser.add_argument('ctf', help='Input ctf.pkl')
    #parser.add_argument('--pose', help='Input pose.pkl') # TODO
    parser.add_argument('--ind', help='Selected indices array (.pkl)')
    parser.add_argument('--full-path', action='store_true', help='Write the full path to particles (default: filename only)')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output .star file')
    return parser

def main(args):
    assert args.o.endswith('.star')
    assert args.particles.endswith('.mrcs'), "Only a single particle stack as an .mrcs is currently supported"
    particles = mrc.parse_mrc(args.particles,lazy=True)[0]
    ctf = utils.load_pkl(args.ctf)
    assert ctf.shape[1] == 9, "Incorrect CTF pkl format"
    assert len(particles) == len(ctf), f"{len(particles)} != {len(ctf)}, Number of particles != number of CTF paraameters"
    log('{} particles'.format(len(particles)))

    if args.ind:
        ind = utils.load_pkl(args.ind)
        log(f'Filtering to {len(ind)} particles')
        ctf = ctf[ind]
    else:
        ind = np.arange(len(particles))

    # _rlnImageName
    ind += 1 # CHANGE TO 1-BASED INDEXING
    image_name = os.path.basename(args.particles) if not args.full_path else args.particles
    names = [f'{i}@{image_name}' for i in ind]
    
    ctf = ctf[:,2:]

    data = {HEADERS[0]:names}
    for i in range(7):
        data[HEADERS[i+1]] = ctf[:,i]
    df = pd.DataFrame(data=data)

    s = starfile.Starfile(HEADERS,df)
    s.write(args.o)

if __name__ == '__main__':
    main(parse_args().parse_args())
