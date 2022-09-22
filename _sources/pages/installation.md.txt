# cryoDRGN Installation

We provide installation instructions assuming an **Anaconda** environment for managing dependencies. Anaconda is a
python package/environment manager which can handle complex dependencies between Python packages through the creation
of python environments. We recommended creating a separate environment for cryodrgn to prevent any conflicts among
dependencies with other software packages.

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

## 1) Install anaconda

- For most platforms, the installation typically consists of a shell script
  (e.g. [Anaconda Installers](https://www.anaconda.com/products/distribution)) that you execute on the command line
  which will prompt you to install and choose a base directory where all the downloaded software and environments will
  go.

  See the official Anaconda documentation and follow their installation instructions
  [here](https://docs.anaconda.com/anaconda/install/).

- Once your anaconda environment is activated, your anaconda environment should be indicated on the command line, e.g.:
    - `(base) $`

## 2) Setting up the cryoDRGN environment

- First, create a new conda environment named `cryodrgn` (or renamed as appropriate):

    ```bash
    (base) $ conda create --name cryodrgn python=3.9
    ```

- Activate the environment. Your command prompt will usually indicate the environment you are in with
  `(environment name)` before the prompt:

    ```bash
    (base) $ conda activate cryodrgn
    (cryodrgn) $
    ```

- Install pytorch and cudatoolkit into your new cryodrgn environment:

    ```bash
    (cryodrgn) $ conda install pytorch cudatoolkit=11.7 -c pytorch
    ```

- Replace the cudatoolkit version with the appropriate version of CUDA installed with the GPU drivers. You can
  check the CUDA version with `nvidia-smi`.
  ```
  +-----------------------------------------------------------------------------+
  | NVIDIA-SMI 515.65.01    Driver Version: 515.65.01    CUDA Version: 11.7     |
  |-------------------------------+----------------------+----------------------+
  | GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
  | Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
  |                               |                      |               MIG M. |
  |===============================+======================+======================|
  |   0  NVIDIA GeForce ...  Off  | 00000000:01:00.0 Off |                  N/A |
  | N/A   41C    P0    N/A /  N/A |      5MiB /  4096MiB |      0%      Default |
  |                               |                      |                  N/A |
  +-------------------------------+----------------------+----------------------+

  +-----------------------------------------------------------------------------+
  | Processes:                                                                  |
  |  GPU   GI   CI        PID   Type   Process name                  GPU Memory |
  |        ID   ID                                                   Usage      |
  |=============================================================================|
  |    0   N/A  N/A      1420      G   /usr/lib/xorg/Xorg                  4MiB |
  +-----------------------------------------------------------------------------+
    ```
  - Don't forget to include `-c pytorch` to get the software from the official pytorch channel
  - To customize the installation line depending on your situation, look at Pytorch's
    [Start locally](https://pytorch.org/get-started/locally/).

- Optional step for older pytorch versions (pre-version 1.6). (Newer versions of pytorch natively support mixed
  precision training.)
    - For accelerated training speeds (~3x faster) on GPUs with tensor cores (Nvidia Volta, Turing, and Ampere GPU
      architectures), install the apex package into the active conda environment. See their official documentation and
      installation instructions [here](https://github.com/NVIDIA/apex#quick-start).

    ```bash
    git clone https://github.com/NVIDIA/apex
    cd apex
    pip install -v --disable-pip-version-check --no-cache-dir ./
    ```

## 3) Install cryoDRGN

- Obtain cryodrgn source code by cloning the git repository, and then doing a `pip install .` in the checkout folder.
  This will also install dependencies that `cryoDRGN` depends on.

```bash
# Clone source code and install
(cryodrgn) $ git clone https://github.com/zhonge/cryodrgn.git
(cryodrgn) $ cd cryodrgn
(cryodrgn) $ pip install .
```

- Alternatively, if you prefer not to use git, you can directly download a ZIP file of the latest release from
  [https://github.com/zhonge/cryodrgn/releases](https://github.com/zhonge/cryodrgn/releases). For example, if the latest
  release is `cryodrgn-1.1.0.zip`:

```bash
(cryodrgn) $ unzip cryodrgn-1.1.0.zip
(cryodrgn) $ cd cryodrgn-1.1.0
(cryodrgn) $ pip install .
```

## 4) Testing the Installation

Once installed, you should be able to call the `cryodrgn` executable and see a list of commands:

```bash
(cryodrgn) $ cryodrgn -h
```

There is a small testing dataset in the source code that you can use to run cryodrgn and verify that all the
dependencies were installed correctly:

```bash
(cryodrgn) $ cd [sourcecode directory]/testing
(cryodrgn) $ ./quicktest.sh
```

It should take ~20 seconds to run and reach a final loss around 0.08 in version 1.0 and 0.03 in version 1.1+. The
output should look something like:

```
+ cryodrgn train_vae data/hand.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/hand_rot.pkl --zdim 10 --pe-type gaussian
2022-09-20 16:30:05     /home/vineetb/.conda/envs/cryodrgn/bin/cryodrgn train_vae data/hand.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/hand_rot.pkl --zdim 10 --pe-type gaussian
2022-09-20 16:30:05     Namespace(particles='/home/vineetb/cryodrgn/cryodrgn/testing/data/hand.mrcs', outdir='/home/vineetb/cryodrgn/cryodrgn/testing/output/toy_recon_vae', zdim=10, poses='/home/vineetb/cryodrgn/cryodrgn/testing/data/hand_rot.pkl', ctf=None, load=None, checkpoint=1, log_interval=1000, verbose=False, seed=0, ind=None, invert_data=True, window=True, window_r=0.85, datadir=None, lazy=False, preprocessed=False, max_threads=16, tilt=None, tilt_deg=45, num_epochs=20, batch_size=8, wd=0, lr=0.0001, beta=None, beta_control=None, norm=None, amp=True, multigpu=False, do_pose_sgd=False, pretrain=1, emb_type='quat', pose_lr=0.0003, qlayers=3, qdim=1024, encode_mode='resid', enc_mask=None, use_real=False, players=3, pdim=1024, pe_type='gaussian', feat_sigma=0.5, pe_dim=None, domain='fourier', activation='relu', func=<function main at 0x7f47d8676790>)
2022-09-20 16:30:06     Use cuda True
2022-09-20 16:30:06     Loading dataset from /home/vineetb/cryodrgn/cryodrgn/testing/data/hand.mrcs
2022-09-20 16:30:06     Loaded 100 64x64 images
2022-09-20 16:30:06     Windowing images with radius 0.85
2022-09-20 16:30:06     Computing FFT
2022-09-20 16:30:06     Spawning 16 processes
2022-09-20 16:30:06     Symmetrizing image data
2022-09-20 16:30:06     Normalized HT by 0 +/- 94.426513671875
2022-09-20 16:30:06     WARNING: No translations provided
2022-09-20 16:30:07     Using circular lattice with radius 32
2022-09-20 16:30:07     HetOnlyVAE(
  (encoder): ResidLinearMLP(
    (main): Sequential(
      (0): Linear(in_features=3208, out_features=1024, bias=True)
      (1): ReLU()
      (2): ResidLinear(
        (linear): Linear(in_features=1024, out_features=1024, bias=True)
      )
      (3): ReLU()
      (4): ResidLinear(
        (linear): Linear(in_features=1024, out_features=1024, bias=True)
      )
      (5): ReLU()
      (6): ResidLinear(
        (linear): Linear(in_features=1024, out_features=1024, bias=True)
      )
      (7): ReLU()
      (8): Linear(in_features=1024, out_features=20, bias=True)
    )
  )
  (decoder): FTPositionalDecoder(
    (decoder): ResidLinearMLP(
      (main): Sequential(
        (0): Linear(in_features=202, out_features=1024, bias=True)
        (1): ReLU()
        (2): ResidLinear(
          (linear): Linear(in_features=1024, out_features=1024, bias=True)
        )
        (3): ReLU()
        (4): ResidLinear(
          (linear): Linear(in_features=1024, out_features=1024, bias=True)
        )
        (5): ReLU()
        (6): ResidLinear(
          (linear): Linear(in_features=1024, out_features=1024, bias=True)
        )
        (7): ReLU()
        (8): Linear(in_features=1024, out_features=2, bias=True)
      )
    )
  )
)
2022-09-20 16:30:07     9814038 parameters in model
2022-09-20 16:30:07     6455316 parameters in encoder
2022-09-20 16:30:07     3358722 parameters in decoder
2022-09-20 16:30:07     Warning: z dimension is not a multiple of 8 -- AMP training speedup is not optimized
2022-09-20 16:30:08     # =====> Epoch: 1 Average gen loss = 1.10873, KLD = 3.386924, total loss = 1.108840; Finished in 0:00:01.346893
2022-09-20 16:30:09     # =====> Epoch: 2 Average gen loss = 0.740441, KLD = 8.101010, total loss = 0.740694; Finished in 0:00:00.542031
2022-09-20 16:30:09     # =====> Epoch: 3 Average gen loss = 0.535575, KLD = 10.675920, total loss = 0.535908; Finished in 0:00:00.541637
2022-09-20 16:30:10     # =====> Epoch: 4 Average gen loss = 0.368407, KLD = 13.592397, total loss = 0.368831; Finished in 0:00:00.541102
2022-09-20 16:30:11     # =====> Epoch: 5 Average gen loss = 0.234344, KLD = 16.974737, total loss = 0.234873; Finished in 0:00:00.544139
2022-09-20 16:30:11     # =====> Epoch: 6 Average gen loss = 0.140454, KLD = 19.307134, total loss = 0.141056; Finished in 0:00:00.541751
2022-09-20 16:30:12     # =====> Epoch: 7 Average gen loss = 0.0899037, KLD = 20.093284, total loss = 0.090530; Finished in 0:00:00.544378
2022-09-20 16:30:13     # =====> Epoch: 8 Average gen loss = 0.0641149, KLD = 20.714445, total loss = 0.064761; Finished in 0:00:00.543761
2022-09-20 16:30:13     # =====> Epoch: 9 Average gen loss = 0.0496119, KLD = 20.798030, total loss = 0.050260; Finished in 0:00:00.545478
2022-09-20 16:30:14     # =====> Epoch: 10 Average gen loss = 0.0408589, KLD = 21.089181, total loss = 0.041516; Finished in 0:00:00.548600
2022-09-20 16:30:15     # =====> Epoch: 11 Average gen loss = 0.0338546, KLD = 21.053582, total loss = 0.034511; Finished in 0:00:00.560902
2022-09-20 16:30:15     # =====> Epoch: 12 Average gen loss = 0.0290218, KLD = 21.509603, total loss = 0.029692; Finished in 0:00:00.549171
2022-09-20 16:30:16     # =====> Epoch: 13 Average gen loss = 0.0252569, KLD = 21.402734, total loss = 0.025924; Finished in 0:00:00.554416
2022-09-20 16:30:17     # =====> Epoch: 14 Average gen loss = 0.0222708, KLD = 21.686829, total loss = 0.022947; Finished in 0:00:00.549451
2022-09-20 16:30:17     # =====> Epoch: 15 Average gen loss = 0.0196031, KLD = 21.829715, total loss = 0.020284; Finished in 0:00:00.547152
2022-09-20 16:30:18     # =====> Epoch: 16 Average gen loss = 0.0175662, KLD = 21.648027, total loss = 0.018241; Finished in 0:00:00.548684
2022-09-20 16:30:19     # =====> Epoch: 17 Average gen loss = 0.0159719, KLD = 21.876881, total loss = 0.016654; Finished in 0:00:00.540562
2022-09-20 16:30:19     # =====> Epoch: 18 Average gen loss = 0.0147737, KLD = 21.754937, total loss = 0.015452; Finished in 0:00:00.541675
2022-09-20 16:30:20     # =====> Epoch: 19 Average gen loss = 0.0133148, KLD = 21.684366, total loss = 0.013991; Finished in 0:00:00.543656
2022-09-20 16:30:20     # =====> Epoch: 20 Average gen loss = 0.0124398, KLD = 21.621814, total loss = 0.013114; Finished in 0:00:00.543454
2022-09-20 16:30:21     Finished in 0:00:15.839427 (0:00:00.791971 per epoch)
```

- You will want to verify that the output contains  `Use cuda True` in the first few lines to ensure that `cryoDRGN`
  will be using your GPU for training.

## Updating cryoDRGN

To update to a later version, you need to obtain the updated software either with `git checkout <version>` or direct
download from [https://github.com/zhonge/cryodrgn](https://github.com/zhonge/cryodrgn), then rerun `$ pip install .`
in your cryodrgn anaconda environment:

```bash
(cryodrgn) $ cd /path/to/repo
(cryodrgn) $ git checkout 1.1.0  # or `git pull origin master` to get the latest sw
(cryodrgn) $ pip install .
```

To keep multiple versions of cryoDRGN in parallel, you will need to create a new anaconda environment and re-install
all the dependencies.

## Known Issues

1. [Jupyter notebook widget not showing up](https://github.com/zhonge/cryodrgn/issues/34)
2. In the jupyter notebook: no attribute `from_dcm` or `from_matrix`:
    - e.g. `AttributeError: type object 'scipy.spatial.transform.rotation.Rotation' has no attribute 'from_dcm'`
    - Make sure you are using scipy version 1.4.0 or later.
    - If you are using scipy version 1.6.0 or later, make sure you are using cryodrgn version 0.3.2 or later.
    - [https://github.com/zhonge/cryodrgn/issues/39](https://github.com/zhonge/cryodrgn/issues/39)
3. In `cryodrgn analyze`: Running UMAP hangs for certain versions of umap
    - Under active investigation [here](https://github.com/zhonge/cryodrgn/issues/53)

If you run into any issues getting cryoDRGN installed, please file a
[github issue](http://www.github.com/zhonge/cryodrgn), including all the commands you used and their output.
