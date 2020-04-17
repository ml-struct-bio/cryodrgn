'''Parse CTF parameters from a cryoSPARC particles.cs file'''

import argparse
import numpy as np
import sys, os
import pickle

from cryodrgn import utils
from cryodrgn import ctf
log = utils.log 

def add_args(parser):
    parser.add_argument('cs', help='Input cryosparc particles.cs file')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output pkl of CTF parameters')
    parser.add_argument('--png', metavar='PNG', type=os.path.abspath, help='Optionally plot the CTF')
    return parser

def main(args):
    assert args.cs.endswith('.cs'), "Input file must be a .cs file"
    assert args.o.endswith('.pkl'), "Output CTF parameters must be .pkl file"

    metadata = np.load(args.cs)
    N = len(metadata)
    log('{} particles'.format(N))

    ctf_params = np.zeros((N, 9))
    ctf_params[:,0] = metadata['blob/shape'][0][0]
    fields = ('blob/psize_A','ctf/df1_A','ctf/df2_A','ctf/df_angle_rad','ctf/accel_kv','ctf/cs_mm','ctf/amp_contrast','ctf/phase_shift_rad')
    for i,f in enumerate(fields):
        ctf_params[:,i+1] = metadata[f]
        if f in ('ctf/df_angle_rad', 'ctf/phase_shift_rad'): # convert to degrees
            ctf_params[:,i+1] *= (180/np.pi) 

    ctf.print_ctf_params(ctf_params[0])
    log('Saving {}'.format(args.o))
    with open(args.o,'wb') as f:
        pickle.dump(ctf_params.astype(np.float32), f)
    if args.png:
        import matplotlib.pyplot as plt
        ctf.plot_ctf(int(ctf_params[0,0]), ctf_params[0,1], ctf_params[0,2:])
        plt.savefig(args.png)
        log(args.png)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    main(add_args(parser).parse_args())
