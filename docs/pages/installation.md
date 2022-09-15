# cryoDRGN installation with anaconda

We provide installation instructions assuming an **Anaconda** environment for managing dependencies. Anaconda is a python package/environment manager which can handle complex dependencies between Python packages through the creation of [environments](http://environments.It). We recommended creating a separate environment for cryodrgn to prevent any conflicts among dependencies with other software packages.

**Compute/hardware requirements:**

- High performance linux workstation or cluster
- NVIDIA GPUs

**Dependencies:**

- python
- pytorch
- cudatoolkit
- numpy
- pandas

**Additional dependencies for visualization:**

- matplotlib
- seaborn
- scipy 1.4.0+
- scikit-learn
- umap
- jupterlab
- ipywidgets
- plotly and cufflinks

The software has been tested on Python 3.7-3.9 and pytorch 1.0-1.7, 1.12.

---

## **Step 1) Install anaconda**

- For linux, the installation typically consists of a shell script (e.g. [`Anaconda3-2019.10-Linux-ppc64le.sh`](http://anaconda3-2019.10-linux-ppc64le.sh/)) that you execute on the command line which will prompt you to install and choose a base directory where all the downloaded software and environments will go.
    - See the official Anaconda documentation and follow their installation instructions [here](https://docs.anaconda.com/anaconda/install/linux/).
- Once your anaconda environment is activated, your anaconda environment should be indicated on the command line, e.g.:
    - `(base) [Tue Feb 02 13:20 zhonge] $`

## **Step 2) Setting up the cryoDRGN environment**

- First, create a new conda environment named `cryodrgn` (or renamed as appropriate):
    
    ```bash
    (base) $ conda create --name cryodrgn python=3.7
    ```
    
- Activate the environment. Your command prompt will usually indicate the environment you are in with `(environment name)` before the prompt:
    
    ```bash
    (base) $ conda activate cryodrgn
    (cryodrgn) $
    ```
    
- Install pytorch and cudatoolkit into your new cryodrgn environment:
    
    ```bash
    (cryodrgn) $ conda install pytorch cudatoolkit=10.2 -c pytorch
    ```
    
    - Replace the cudatoolkit version with the appropriate version of CUDA installed with the GPU drivers (you can check the CUDA version with `nvidia-smi` , [example here](https://varhowto.com/check-cuda-version/#Method_2_%25E2%2580%2594_Check_CUDA_version_by_nvidiasmi_from_NVIDIA_Linux_driver))
    - Don't forget to include `-c pytorch` to get the software from the official pytorch channel
    - For more detailed installation instructions, see the official pytorch documentation [here](https://pytorch.org/get-started/locally/).
- Install other cryodrgn dependencies. These are common python packages, and anaconda will automatically find the appropriate version that is compatible with your system.
    
    ```python
    (cryodrgn) $ conda install pandas
    (cryodrgn) $ conda install seaborn scikit-learn 
    (cryodrgn) $ conda install umap-learn jupyterlab ipywidgets cufflinks-py "nodejs>=15.12.0" -c conda-forge
    (cryodrgn) $ jupyter labextension install @jupyter-widgets/jupyterlab-manager --no-build
    (cryodrgn) $ jupyter labextension install jupyterlab-plotly --no-build
    (cryodrgn) $ jupyter labextension install plotlywidget --no-build
    (cryodrgn) $ jupyter lab build
    ```
    
- Optional step for older pytorch versions (pre-version 1.6). (Newer versions of pytorch natively support mixed precision training.)
    - For accelerated training speeds (~3x faster) on GPUs with tensor cores (Nvidia Volta, Turing, and Ampere GPU architectures), install the apex package into the conda environment
        - See their official documentation and installation instructions here: [https://github.com/NVIDIA/apex#quick-start](https://github.com/NVIDIA/apex#quick-start)
    
    ```python
    git clone https://github.com/NVIDIA/apex
    cd apex
    pip install -v --disable-pip-version-check --no-cache-dir ./
    ```
    

## Step 3) Install cryoDRGN

- Obtain cryodrgn source code by cloning the git repository (a github account is required):

```python
# Clone source code and install
(cryodrgn) $ git clone https://github.com/zhonge/cryodrgn.git
(cryodrgn) $ cd cryodrgn
(cryodrgn) $ git checkout 1.1.09 # Or latest released version
(cryodrgn) $ python setup.py install
```

- Alternatively, if you prefer not to use git, you can directly download a ZIP file of the source code from a browser [https://github.com/zhonge/cryodrgn/releases](https://github.com/zhonge/cryodrgn/releases).

```python
(cryodrgn) $ unzip cryodrgn-1.1.0.zip
(cryodrgn) $ cd cryodrgn-1.1.0
(cryodrgn) $ python setup.py install
```

## Step 4) Testing the Installation

Once installed, you should be able to call the `cryodrgn` executable and see a list of commands:

```python
(cryodrgn) $ cryodrgn -h
```

There is a small testing dataset in the source code that you can use to run cryodrgn and verify that all the dependencies were installed correctly: 

```python
(cryodrgn) $ cd [sourcecode directory]/testing
(cryodrgn) $ ./quicktest.sh
```

It should take ~20 seconds to run and reach a final loss around 0.08 in version 1.0 and 0.03 in version 1.1+. The output should look something like:

- `[quicktest.sh](http://quicktest.sh)` output
    
    ```python
    ++ cryodrgn train_vae data/hand.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/hand_rot.pkl --zdim 10
    2021-02-02 14:04:10     /nobackup/users/zhonge/anaconda3/envs/cryodrgn4/bin/cryodrgn train_vae data/hand.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/hand_rot.pkl --zdim 10
    2021-02-02 14:04:10     Namespace(activation='relu', amp=False, batch_size=8, beta=None, beta_control=None, checkpoint=1, ctf=None, datadir=None, do_pose_sgd=False, domain='fourier', emb_type='quat', enc_mask=None, encode_mode='resid', func=<function main at 0x2000ab053840>, ind=None, invert_data=True, lazy=False, load=None, log_interval=1000, lr=0.0001, multigpu=False, norm=None, num_epochs=20, outdir='/nobackup/users/zhonge/dev/cryodrgn/master/testing/output/toy_recon_vae', particles='/nobackup/users/zhonge/dev/cryodrgn/master/testing/data/hand.mrcs', pdim=256, pe_dim=None, pe_type='geom_lowf', players=3, pose_lr=0.0003, poses='/nobackup/users/zhonge/dev/cryodrgn/master/testing/data/hand_rot.pkl', pretrain=1, qdim=256, qlayers=3, relion31=False, seed=0, tilt=None, tilt_deg=45, use_real=False, verbose=False, wd=0, window=True, zdim=10)
    2021-02-02 14:04:10     Use cuda True
    2021-02-02 14:04:10     Loaded 100 64x64 images
    2021-02-02 14:04:10     Normalized HT by 0 +/- 94.426513671875
    2021-02-02 14:04:10     WARNING: No translations provided
    2021-02-02 14:04:22     Using circular lattice with radius 32
    2021-02-02 14:04:22     HetOnlyVAE(
      (encoder): ResidLinearMLP(
        (main): Sequential(
          (0): Linear(in_features=3208, out_features=256, bias=True)
          (1): ReLU()
          (2): ResidLinear(
            (linear): Linear(in_features=256, out_features=256, bias=True)
          )
          (3): ReLU()
          (4): ResidLinear(
            (linear): Linear(in_features=256, out_features=256, bias=True)
          )
          (5): ReLU()
          (6): ResidLinear(
            (linear): Linear(in_features=256, out_features=256, bias=True)
          )
          (7): ReLU()
          (8): Linear(in_features=256, out_features=20, bias=True)
        )
      )
      (decoder): FTPositionalDecoder(
        (decoder): ResidLinearMLP(
          (main): Sequential(
            (0): Linear(in_features=202, out_features=256, bias=True)
            (1): ReLU()
            (2): ResidLinear(
              (linear): Linear(in_features=256, out_features=256, bias=True)
            )
            (3): ReLU()
            (4): ResidLinear(
              (linear): Linear(in_features=256, out_features=256, bias=True)
            )
            (5): ReLU()
            (6): ResidLinear(
              (linear): Linear(in_features=256, out_features=256, bias=True)
            )
            (7): ReLU()
            (8): Linear(in_features=256, out_features=2, bias=True)
          )
        )
      )
    )
    2021-02-02 14:04:22     1273878 parameters in model
    2021-02-02 14:04:26     # =====> Epoch: 1 Average gen loss = 1.20397, KLD = 1.443059, total loss = 1.204012; Finished in 0:00:03.281083
    2021-02-02 14:04:26     # =====> Epoch: 2 Average gen loss = 1.10656, KLD = 3.117473, total loss = 1.106653; Finished in 0:00:00.199297
    2021-02-02 14:04:26     # =====> Epoch: 3 Average gen loss = 1.03491, KLD = 4.985058, total loss = 1.035065; Finished in 0:00:00.198787
    2021-02-02 14:04:27     # =====> Epoch: 4 Average gen loss = 0.954972, KLD = 6.849713, total loss = 0.955186; Finished in 0:00:00.198769
    2021-02-02 14:04:27     # =====> Epoch: 5 Average gen loss = 0.871234, KLD = 11.276586, total loss = 0.871585; Finished in 0:00:00.199222
    2021-02-02 14:04:27     # =====> Epoch: 6 Average gen loss = 0.784677, KLD = 14.761210, total loss = 0.785138; Finished in 0:00:00.198983
    2021-02-02 14:04:28     # =====> Epoch: 7 Average gen loss = 0.694899, KLD = 18.206867, total loss = 0.695466; Finished in 0:00:00.199065
    2021-02-02 14:04:28     # =====> Epoch: 8 Average gen loss = 0.60537, KLD = 20.765765, total loss = 0.606018; Finished in 0:00:00.199060
    2021-02-02 14:04:28     # =====> Epoch: 9 Average gen loss = 0.522229, KLD = 23.063340, total loss = 0.522948; Finished in 0:00:00.200438
    2021-02-02 14:04:29     # =====> Epoch: 10 Average gen loss = 0.442045, KLD = 24.941707, total loss = 0.442823; Finished in 0:00:00.198947
    2021-02-02 14:04:29     # =====> Epoch: 11 Average gen loss = 0.37073, KLD = 25.814933, total loss = 0.371535; Finished in 0:00:00.199022
    2021-02-02 14:04:29     # =====> Epoch: 12 Average gen loss = 0.306593, KLD = 27.169124, total loss = 0.307440; Finished in 0:00:00.199095
    2021-02-02 14:04:30     # =====> Epoch: 13 Average gen loss = 0.254798, KLD = 27.483357, total loss = 0.255655; Finished in 0:00:00.200581
    2021-02-02 14:04:30     # =====> Epoch: 14 Average gen loss = 0.20721, KLD = 28.741093, total loss = 0.208106; Finished in 0:00:00.198856
    2021-02-02 14:04:31     # =====> Epoch: 15 Average gen loss = 0.16829, KLD = 28.868599, total loss = 0.169189; Finished in 0:00:00.202039
    2021-02-02 14:04:31     # =====> Epoch: 16 Average gen loss = 0.138856, KLD = 29.797036, total loss = 0.139785; Finished in 0:00:00.199184
    2021-02-02 14:04:31     # =====> Epoch: 17 Average gen loss = 0.117319, KLD = 30.307774, total loss = 0.118264; Finished in 0:00:00.198926
    2021-02-02 14:04:32     # =====> Epoch: 18 Average gen loss = 0.100041, KLD = 31.086511, total loss = 0.101010; Finished in 0:00:00.199008
    2021-02-02 14:04:32     # =====> Epoch: 19 Average gen loss = 0.0883892, KLD = 31.568129, total loss = 0.089373; Finished in 0:00:00.200003
    2021-02-02 14:04:32     # =====> Epoch: 20 Average gen loss = 0.079563, KLD = 31.856020, total loss = 0.080556; Finished in 0:00:00.199169
    2021-02-02 14:04:33     Finsihed in 0:00:22.683726 (0:00:01.134186 per epoch)
    ```
    
- Note that the output should contain  `Use cuda True` in the first few lines

## Updating cryoDRGN versions

To update to a later version, you need to obtain the updated software either with `git checkout VERSION` or direct download from [https://github.com/zhonge/cryodrgn](https://github.com/zhonge/cryodrgn), then rerun `$ python [setup.py](http://setup.py) install` in your cryodrgn anaconda environment:

```python
(cryodrgn) $ cd /path/to/repo
(cryodrgn) $ git checkout 1.1.0 # or `git pull origin master` to get the latest sw
(cryodrgn) $ python setup.py install
```

To keep multiple versions of cryoDRGN in parallel, you will need to create a new anaconda environment and re-install all the dependencies.

## Known Issues

1. [Jupyter notebook widget not showing up](https://github.com/zhonge/cryodrgn/issues/34)
2. In the jupyter notebook: no attribute `from_dcm` or `from_matrix`:
    - e.g. `AttributeError: type object 'scipy.spatial.transform.rotation.Rotation' has no attribute 'from_dcm'`
    - Make sure you are using scipy version 1.4.0 or later.
    - If you are using scipy version 1.6.0 or later, make sure you are using cryodrgn version 0.3.2 or later.
    - [https://github.com/zhonge/cryodrgn/issues/39](https://github.com/zhonge/cryodrgn/issues/39)
3. In `cryodrgn analyze`: Running UMAP hangs for certain versions of umap
    - Under active investigation [here](https://github.com/zhonge/cryodrgn/issues/53)

If you run into any issues getting cryoDRGN installed, please file a [github issue](http://www.github.com/zhonge/cryodrgn), including all the commands you used and their output.