#!/usr/bin/env python

from distutils.core import setup

setup(name='cryodrgn',
      version='0.2.0',
      description='cryoDRGN heterogeneous reconstruction',
      author='Ellen Zhong',
      author_email='zhonge@mit.edu',
      url='https://github.com/zhonge/cryodrgn',
      license='GPLv3',
      packages=['cryodrgn'],
      install_requires=[
        'torch>=1.0.0',
        'torchvision',
        'pandas',
        ]
     )
