'''Parse 3D alignments from cryosparc .cs metafiles'''

import argparse
import numpy as np
import sys, os
import pickle
sys.path.insert(0,'{}/../lib-python'.format(os.path.dirname(os.path.abspath(__file__))))
import utils
log = utils.log

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Cryosparc .cs file')
    parser.add_argument('-o', help='Output prefix for appending .rot.pkl and .trans.pkl')
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
    body = [x for x in body if x.strip()] # remove any empty lines
    body = [x.split() for x in body]
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

def extract(header, body, keys):
    ind = [header.index(x) for x in keys]
    euler = [[x[i] for i in ind] for x in body]
    euler = np.asarray(euler, dtype=np.float32)
    return euler

def main(args):
    header, body = parse_star(args.input)
    log('{} rows'.format(len(body)))
    
    # parse rotations
    keys = ('_rlnAngleRot','_rlnAngleTilt','_rlnAnglePsi')
    euler = extract(header, body, keys)
    log(euler[0])
    rot = np.asarray([utils.R_from_relion(*x) for x in euler])
    log(rot[0])
    log(rot.shape)

    # parse translations
    keys = ('_rlnOriginX','_rlnOriginY')
    trans = extract(header, body, keys)
    log(trans[0])
    log(trans.shape)
    
    # write output
    out_rot = '{}.rot.pkl'.format(args.o)
    log('Writing {}'.format(out_rot))
    with open(out_rot,'wb') as f:
        pickle.dump(rot,f)
    out_trans = '{}.trans.pkl'.format(args.o)
    log('Writing {}'.format(out_trans))
    with open(out_trans,'wb') as f:
        pickle.dump(trans,f)

if __name__ == '__main__':
    main(parse_args().parse_args())
