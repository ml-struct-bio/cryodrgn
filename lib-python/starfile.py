'''
Lightweight parser for starfiles
'''

import numpy as np
import pandas as pd
from datetime import datetime as dt

class Starfile():
    
    def __init__(self, headers, df):
        self.headers = headers
        self.df = df

    @classmethod
    def load_starfile(self, starfile):
        f = open(starfile,'r')
        # get to data block
        while 1:
            for line in f:
                if line.startswith('data_'):
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
        # assume the rest is the body
        headers = [h.strip().split()[0] for h in headers]
        body = [line] + f.readlines()
        # remove last line of body if empty
        if body[-1].strip() == '':
            body = body[:-1]
        # put data into an array and instantiate as dataframe
        words = [l.strip().split() for l in body]
        words = np.array(words)
        data = {h:words[:,i] for i,h in enumerate(headers)}
        df = pd.DataFrame(data=d)
        return self(headers, df)

    def write(self, outstar):
        f = open(outstar,'w')
        f.write('# Created {}\n'.format(dt.now()))
        f.write('\n')
        f.write('data_\n\n')
        f.write('loop_\n')
        f.write('\n'.join(self.headers))
        f.write('\n')
        for i in range(len(self.df)):
            f.write(' '.join(self.df.loc[i]))
            f.write('\n')
        #f.write('\n'.join([' '.join(self.df.loc[i]) for i in range(len(self.df))]))
