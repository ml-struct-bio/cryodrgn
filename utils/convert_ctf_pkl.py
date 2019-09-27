'''Parse CTF params from a RELION starfile'''

import argparse
import numpy as np
import sys, os
import pickle

sys.path.insert(0, '{}/../lib-python'.format(os.path.dirname(os.path.realpath(__file__))))
import utils
import ctf
log = utils.log 

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('star', help='Input')
    parser.add_argument('-N', type=int, required=True, help='Number of particles in image stack')
    parser.add_argument('--Apix', type=float, required=True, help='Angstroms per pixel')
    parser.add_argument('-D', type=int, help='Image size (width of image in pixels)')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output pkl with CTF params')
    parser.add_argument('--png', type=os.path.abspath)
    return parser

def parse_star(starfile):
    f = open(starfile,'r')
    for line in f:
        if line.startswith('loop_'):
            break
    lines = f.readlines()
    header = []
    i = 0
    for l in lines:
        if l.startswith('_rln'):
            header.append(l.strip().split()[0])
            i += 1
        else:
            break
    body = lines[i:]
    body = [x for x in body if x != '\n'] # remove any empty lines
    return header, body

def parse_ctf(starfile, N):
    header, body = parse_star(starfile)
    assert len(body) == N
    body = [x.split() for x in body]
    ind = (header.index(x) for x in ('_rlnDefocusU', '_rlnDefocusV', '_rlnDefocusAngle', '_rlnVoltage', '_rlnSphericalAberration', '_rlnAmplitudeContrast'))
    ind = tuple(ind)
    ctf_params = [[x[i] for i in ind] for x in body]
    ctf_params = np.asarray(ctf_params)
    return ctf_params

def main(args):
    ctf_params = np.zeros((args.N, 7))
    ctf_params[:,0] = args.Apix
    ctf_params[:,1:] = parse_ctf(args.star, args.N)
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
