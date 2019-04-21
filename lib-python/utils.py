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

def load_pkl(pkl):
    with open(pkl,'rb') as f:
        x = pickle.load(f)
    return x

def load_angles(pkl):
    return load_pkl(pkl)

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

