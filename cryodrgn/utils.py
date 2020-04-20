from datetime import datetime as dt
import os, sys
import numpy as np
import pickle

_verbose = False

def log(msg):
    print('{}     {}'.format(dt.now().strftime('%Y-%m-%d %H:%M:%S'), msg))
    sys.stdout.flush()

def vlog(msg):
    if _verbose:
        print('{}     {}'.format(dt.now().strftime('%Y-%m-%d %H:%M:%S'), msg))
        sys.stdout.flush()

def flog(msg, outfile):
    msg = '{}     {}'.format(dt.now().strftime('%Y-%m-%d %H:%M:%S'), msg)
    print(msg)
    sys.stdout.flush()
    try:
        with open(outfile,'a') as f:
            f.write(msg+'\n')
    except Exception as e:
        log(e)

def load_pkl(pkl):
    with open(pkl,'rb') as f:
        x = pickle.load(f)
    return x

def save_pkl(data, out_pkl, append='False'):
    mode = 'wb' if append is False else 'ab'
    if mode == 'wb' and os.path.exists(out_pkl):
        vlog(f'Warning: {out_pkl} already exists. Overwriting.')
    with open(out_pkl, mode) as f:
        x = pickle.dump(data, f)

def R_from_eman(a,b,y):
    a *= np.pi/180.
    b *= np.pi/180.
    y *= np.pi/180.
    ca, sa = np.cos(a), np.sin(a)
    cb, sb = np.cos(b), np.sin(b)
    cy, sy = np.cos(y), np.sin(y)
    Ra = np.array([[ca,-sa,0],[sa,ca,0],[0,0,1]])
    Rb = np.array([[1,0,0],[0,cb,-sb],[0,sb,cb]])
    Ry = np.array(([cy,-sy,0],[sy,cy,0],[0,0,1]))
    R = np.dot(np.dot(Ry,Rb),Ra)
    # handling EMAN convention mismatch for where the origin of an image is (bottom right vs top right)
    R[0,1] *= -1
    R[1,0] *= -1
    R[1,2] *= -1
    R[2,1] *= -1
    return R

def R_from_relion(a,b,y):
    a *= np.pi/180.
    b *= np.pi/180.
    y *= np.pi/180.
    ca, sa = np.cos(a), np.sin(a)
    cb, sb = np.cos(b), np.sin(b)
    cy, sy = np.cos(y), np.sin(y)
    Ra = np.array([[ca,-sa,0],[sa,ca,0],[0,0,1]])
    Rb = np.array([[cb,0,-sb],[0,1,0],[sb,0,cb]])
    Ry = np.array(([cy,-sy,0],[sy,cy,0],[0,0,1]))
    R = np.dot(np.dot(Ry,Rb),Ra)
    R[0,1] *= -1
    R[1,0] *= -1
    R[1,2] *= -1
    R[2,1] *= -1
    return R

def xrot(tilt_deg):
    '''Return rotation matrix associated with rotation over the x-axis'''
    theta = tilt_deg*np.pi/180
    tilt = np.array([[1.,0.,0.],
                     [0, np.cos(theta), -np.sin(theta)],
                     [0, np.sin(theta), np.cos(theta)]])
    return tilt
    
def zero_sphere(vol):
    '''Zero values of @vol outside the sphere'''
    assert len(set(vol.shape)) == 1, 'volume must be a cube'
    D = vol.shape[0]
    xx = np.linspace(-1, 1, D, endpoint=True if D % 2 == 1 else False)
    z,y,x = np.meshgrid(xx,xx,xx)
    coords = np.stack((x,y,z),-1)
    r = np.sum(coords**2,axis=-1)**.5
    vlog('Zeroing {} pixels'.format(len(np.where(r>1)[0])))
    vol[np.where(r>1)] = 0
    return vol

