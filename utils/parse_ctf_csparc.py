'''Parse CTF parameters from a cryoSPARC particles.cs file'''

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
    parser.add_argument('cs', help='Input cryosparc particles.cs file')
    parser.add_argument('-N', type=int, required=True, help='Number of particles in image stack')
    parser.add_argument('--Apix', type=float, required=True, help='Angstroms per pixel')
    parser.add_argument('-D', type=int, help='Image size in pixels')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output pkl with CTF params')
    parser.add_argument('--png', metavar='PNG', type=os.path.abspath, help='Optionally plot the CTF')
    return parser

#      dtype=[('uid', '<u8'), ('blob/path', 'S89'), ('blob/idx', '<u4'), ('blob/shape', '<u4', (2,)), ('blob/psize_A', '<f4'), ('blob/sign', '<f4'), ('ctf/type', 'S9'), ('ctf/exp_group_id', '<u4'), ('ctf/accel_kv', '<f4'), ('ctf/cs_mm', '<f4'), ('ctf/amp_contrast', '<f4'), ('ctf/df1_A', '<f4'), ('ctf/df2_A', '<f4'), ('ctf/df_angle_rad', '<f4'), ('ctf/phase_shift_rad', '<f4'), ('ctf/scale', '<f4'), ('ctf/scale_const', '<f4'), ('alignments3D/split', '<u4'), ('alignments3D/shift', '<f4', (2,)), ('alignments3D/pose', '<f4', (3,)), ('alignments3D/psize_A', '<f4'), ('alignments3D/error', '<f4'), ('alignments3D/error_min', '<f4'), ('alignments3D/resid_pow', '<f4'), ('alignments3D/slice_pow', '<f4'), ('alignments3D/image_pow', '<f4'), ('alignments3D/cross_cor', '<f4'), ('alignments3D/alpha', '<f4'), ('alignments3D/weight', '<f4'), ('alignments3D/pose_ess', '<f4'), ('alignments3D/shift_ess', '<f4'), ('alignments3D/class_posterior', '<f4'), ('alignments3D/class', '<u4'), ('alignments3D/class_ess', '<f4')])

def main(args):
    assert args.o.endswith('.pkl'), "Output CTF parameters must be .pkl file"
    metadata = np.load(args.cs)
    assert len(metadata) == args.N
    ctf_params = np.zeros((len(metadata), 8))
    fields = ('blob/psize_A','ctf/df1_A','ctf/df2_A','ctf/df_angle_rad','ctf/accel_kv','ctf/cs_mm','ctf/amp_contrast','ctf/phase_shift_rad')
    for i,f in enumerate(fields):
        ctf_params[:,i] = metadata[f]
        if f in ('ctf/df_angle_rad', 'ctf/phase_shift_rad'): # convert to degrees
            ctf_params[:,i] *= (180/np.pi) 

    # check consistency of boxsize
    assert metadata['blob/shape'][0][0]*metadata['blob/psize_A'][0] == args.D*args.Apix
    ctf_params[:,0] = args.Apix
    ctf.print_ctf_params(ctf_params[0])
    log('Saving {}'.format(args.o))
    with open(args.o,'wb') as f:
        pickle.dump(ctf_params.astype(np.float32), f)
    if args.png:
        import matplotlib.pyplot as plt
        import seaborn as sns
        D = args.D
        Apix = ctf_params[0][0]
        freqs = np.stack(np.meshgrid(np.linspace(-.5,.5,D,endpoint=False),np.linspace(-.5,.5,D,endpoint=False)),-1)/Apix
        freqs = freqs.reshape(-1,2)
        c = ctf.compute_ctf_np(freqs, *ctf_params[0,1:])
        sns.heatmap(c.reshape(D, D))
        plt.savefig(args.png)
        log(args.png)
    

if __name__ == '__main__':
    main(parse_args().parse_args())
