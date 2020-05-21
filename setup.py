#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='cryodrgn',
      version='0.2.1',
      description='cryoDRGN heterogeneous reconstruction',
      author='Ellen Zhong',
      author_email='zhonge@mit.edu',
      url='https://github.com/zhonge/cryodrgn',
      license='GPLv3',
      packages=find_packages(),
      entry_points={
          "console_scripts": [
            "cryodrgn = cryodrgn.__main__:main",
            ],
      },
      include_package_data = True,
      install_requires=[
        'torch>=1.0.0',
        'pandas',
        'numpy',
        'matplotlib',
        'scipy>=1.3.1',
        'scikit-learn',
        'seaborn',
        'cufflinks',
        'jupyterlab',
        'umap-learn',
        'ipywidgets'
        ]
     )
