'''Parse CTF parameters from a RELION .star file'''

import argparse
import numpy as np
import sys, os
import pickle

sys.path.insert(0, '{}/../lib-python'.format(os.path.dirname(os.path.realpath(__file__))))
import utils
import starfile
import ctf
log = utils.log 

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('star', help='Input')
    parser.add_argument('-N', type=int, required=True, help='Number of particles in image stack')
    parser.add_argument('--Apix', type=float, required=True, help='Angstroms per pixel')
    parser.add_argument('-D', type=int, help='Image size in pixels')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output pkl with CTF params')
    parser.add_argument('--png', metavar='PNG', type=os.path.abspath, help='Optionally plot the CTF')
    return parser

def main(args):
    assert args.o.endswith('.pkl'), "Output CTF parameters must be .pkl file"
    
    s = starfile.Starfile.load(args.star)
    assert args.N == len(s.df)
    log('{} particles'.format(args.N))

    ctf_params = np.zeros((args.N, 8))
    ctf_params[:,0] = args.Apix
    for i,header in enumerate(['_rlnDefocusU', '_rlnDefocusV', '_rlnDefocusAngle', '_rlnVoltage', '_rlnSphericalAberration', '_rlnAmplitudeContrast', '_rlnPhaseShift']):
        ctf_params[:,i+1] = s.df[header]
    ctf.print_ctf_params(ctf_params[0])
    log('Saving {}'.format(args.o))
    with open(args.o,'wb') as f:
        pickle.dump(ctf_params.astype(np.float32), f)
    if args.png:
        assert args.D, 'Need image size to plot CTF'
        import matplotlib.pyplot as plt
        import seaborn as sns
        freqs = np.stack(np.meshgrid(np.linspace(-.5,.5,args.D,endpoint=False),np.linspace(-.5,.5,args.D,endpoint=False)),-1)/args.Apix
        freqs = freqs.reshape(-1,2)
        c = ctf.compute_ctf_np(freqs, *ctf_params[0,1:])
        sns.heatmap(c.reshape(args.D, args.D))
        plt.savefig(args.png)
        log(args.png)
    

if __name__ == '__main__':
    main(parse_args().parse_args())
