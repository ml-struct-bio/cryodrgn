'''
Lightweight parser for starfiles
'''

import numpy as np
import pandas as pd
from datetime import datetime as dt
import os

from . import mrc
from .mrc import LazyImage

class Starfile():
    
    def __init__(self, headers, df, headers_optics=None, df_optics=None):
        assert headers == list(df.columns), f'{headers} != {df.columns}'
        self.headers = headers
        self.df = df
        if headers_optics is not None: # file is recognized as a RELION 3.1 star file
            self.headers_optics = headers_optics
            self.df_optics = df_optics

    @classmethod
    def load(self, starfile, relion31=False):
        def parse_starblock_to_df(filehandle, BLOCK):
            # get to block
            while 1:
                for line in f:
                    if line.startswith(BLOCK):
                        assert line.strip().split()[0] == BLOCK, f"Expected block header: {BLOCK} but found: {line.strip().split()[0]}. Is this a RELION 3.1 starfile?"
                        break
                break
            # get to header loop
            while 1:
                for line in f:
                    if line.startswith('loop_'):
                        break
                break
            # get list of column headers
            while 1:
                headers = []
                for line in f:
                    if line.startswith('_'):
                        headers.append(line)
                    else:
                        break
                break
            # assume all subsequent lines until empty line is the body
            headers = [h.strip().split()[0] for h in headers]
            assert headers != [], "Could not find headers in file. Please double check if this is a RELION 3.1 starfile."
            if line.strip() != '':
                body = [line]
            else:
                body = []
            for line in f:
                if line.strip() == '':
                    break
                body.append(line)
            # put data into an array and instantiate as dataframe
            words = [l.strip().split() for l in body]
            words = np.array(words)
            assert words.shape[1] == len(headers), f"Error in parsing. Number of columns {words.shape[1]} != number of headers {len(headers)}"
            data = {h: words[:, i] for i, h in enumerate(headers)}
            df = pd.DataFrame(data=data)
            return headers, df

        f = open(starfile,'r')

        # parse data_optics block if relion3.1 starfile
        if relion31:
            BLOCK = 'data_optics'
            headers_optics, df_optics = parse_starblock_to_df(f, BLOCK)

        # parse data block
        BLOCK = 'data_particles' if relion31 else 'data_'
        headers, df = parse_starblock_to_df(f, BLOCK)

        if relion31:
            return self(headers, df, headers_optics, df_optics)
        else:
            return self(headers, df)

    def write(self, outstar, relion31=False):
        f = open(outstar,'w')
        f.write('# Created {}\n'.format(dt.now()))
        f.write('\n')
        if relion31:
            f.write('data_optics\n\n')
            f.write('loop_\n')
            f.write('\n'.join(self.headers_optics))
            f.write('\n')
            for i in self.df_optics.index:
                f.write(' '.join([str(v) for v in self.df_optics.loc[i]]))
                f.write('\n\n')
        f.write('data_particles\n\n') if relion31 else f.write('data_\n\n')
        f.write('loop_\n')
        f.write('\n'.join(self.headers))
        f.write('\n')
        for i in self.df.index:
            f.write(' '.join([str(v) for v in self.df.loc[i]]))
            f.write('\n')
        #f.write('\n'.join([' '.join(self.df.loc[i]) for i in range(len(self.df))]))

    def get_particles(self, datadir=None, lazy=True):
        '''
        Return particles of the starfile

        Input:
            datadir (str): Overwrite base directories of particle .mrcs
                Tries both substituting the base path and prepending to the path
            If lazy=True, returns list of LazyImage instances, else np.array
        '''
        particles = self.df['_rlnImageName']

        # format is index@path_to_mrc
        particles = [x.split('@') for x in particles]
        ind = [int(x[0])-1 for x in particles] # convert to 0-based indexing
        mrcs = [x[1] for x in particles]
        if datadir is not None:
            mrcs = prefix_paths(mrcs, datadir)
        for path in set(mrcs):
            assert os.path.exists(path), f'{path} not found'
        header = mrc.parse_header(mrcs[0])
        D = header.D # image size along one dimension in pixels
        dtype = header.dtype
        stride = dtype().itemsize*D*D
        dataset = [LazyImage(f, (D,D), dtype, 1024+ii*stride) for ii,f in zip(ind, mrcs)]
        if not lazy:
            dataset = np.array([x.get() for x in dataset])
        return dataset

def prefix_paths(mrcs, datadir):
    mrcs1 = ['{}/{}'.format(datadir, os.path.basename(x)) for x in mrcs]
    mrcs2 = ['{}/{}'.format(datadir, x) for x in mrcs]
    try:
        for path in set(mrcs1):
            assert os.path.exists(path)
        mrcs = mrcs1
    except:
        for path in set(mrcs2):
            assert os.path.exists(path), f'{path} not found'
        mrcs = mrcs2
    return mrcs

def csparc_get_particles(csfile, datadir=None, lazy=True):
    metadata = np.load(csfile)
    ind = metadata['blob/idx'] # 0-based indexing
    mrcs = metadata['blob/path'].astype(str).tolist()
    if datadir is not None:
        mrcs = prefix_paths(mrcs, datadir)
    for path in set(mrcs):
        assert os.path.exists(path), f'{path} not found'
    D = metadata[0]['blob/shape'][0]
    dtype = np.float32
    stride = np.float32().itemsize*D*D
    dataset = [LazyImage(f, (D,D), dtype, 1024+ii*stride) for ii,f in zip(ind, mrcs)]
    if not lazy:
        dataset = np.array([x.get() for x in dataset])
    return dataset




