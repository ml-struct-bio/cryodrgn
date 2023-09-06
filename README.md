[![CI](https://github.com/zhonge/cryodrgn/actions/workflows/main.yml/badge.svg)](https://github.com/zhonge/cryodrgn/actions/workflows/main.yml)

# :snowflake::dragon: cryoDRGN: Deep Reconstructing Generative Networks for cryo-EM and cryo-ET heterogeneous reconstruction

CryoDRGN is a neural network based algorithm for heterogeneous cryo-EM reconstruction. In particular, the method models a *continuous* distribution over 3D structures by using a neural network based representation for the volume.

## Documentation:

The latest documentation for cryoDRGN is available [on gitbook](https://ez-lab.gitbook.io/cryodrgn/), including an overview and walkthrough of cryoDRGN installation, training and analysis. A brief quick start is provided below.

For any feedback, questions, or bugs, please file a Github issue, start a Github discussion, or email the [google group](https://groups.google.com/g/cryodrgn).

## New in Version 3.x

The official [cryoDRGN-ET](https://www.biorxiv.org/content/10.1101/2023.08.18.553799v1) release for heterogeneous subtomogram analysis. 

* [NEW] Heterogeneous reconstruction of subtomograms. See documentation [on gitbook](https://ez-lab.gitbook.io/cryodrgn/)
* [NEW] `cryodrgn direct_traversal` for making movies
* Updated `cryodrgn backproject_voxel` for voxel-based homogeneous reconstruction
* Major refactor of dataset loading for handling large datasets

### Previous versions

<details><summary>Version 2.3</summary>
	
* Model configuration files are now saved as human-readable config.yaml files (https://github.com/zhonge/cryodrgn/issues/235)
* Fix machine stamp in output .mrc files for better compatibility with downstream tools (https://github.com/zhonge/cryodrgn/pull/260)
* Better documentation of help flags in ab initio reconstruction tools (https://github.com/zhonge/cryodrgn/issues/258)
* [FIX] By default, window images in `cryodrgn abinit_homo` (now consistent with other reconstruction tools) (https://github.com/zhonge/cryodrgn/issues/258)
* [FIX] Reduce memory usage when using `--preprocessed` and `--ind` (https://github.com/zhonge/cryodrgn/pull/272)

</details>

<details><summary>Version 2.2</summary>

* [NEW] Tools for ab initio homogeneous and heterogeneous reconstruction:

```
(cryodrgn) $ cryodrgn abinit_homo -h
(cryodrgn) $ cryodrgn abinit_het -h
```

* [NEW] Utils function for writing cryoSPARC `.cs`/`.csg` files [to reimport data into cryoSPARC](https://github.com/zhonge/cryodrgn/issues/150#issuecomment-1465094751):

```
(cryodrgn) $ cryodrgn_utils write_cs
```

* [Improved plotting](https://github.com/zhonge/cryodrgn/issues/219) in `cryodrgn analyze`

* Many codebase improvements with open-source software development practices (e.g. continuous integration tests, black, flake8, pyright, logging, and PyPi packaging).

* Note: we are working on a major refactor of data loading for handling large datasets for the next minor version (v2.4). This will entail an API change for the mrc.py library module

</details>

<details><summary>Version 1.1.x</summary>

Updated default parameters for `cryodrgn train_vae` with modified positional encoding, larger model architecture, and accelerated mixed-precision training turned on by default:
* Mixed precision training is now turned on by default (Use `--no-amp` to revert to single precision training)
* Encoder/decoder architecture is now 1024x3 by default (Use `--enc-dim 256` and `--dec-dim 256` to revert)
* Gaussian Fourier featurization for faster training and higher resolution density maps (Use `--pe-type geom_lowf` to revert)

</details>


<details><summary>Version 1.0.x</summary>

The official version 1.0 release. This version introduces several new tools for analysis of the reconstructed ensembles, and adds functionality for calling utility scripts with `cryodrgn_utils <command>`.

* NEW: `cryodrgn analyze_landscape` and `cryodrgn analyze_landscape_full` for automatic assignment of classes and conformational landscape visualization. Documentation for this new feature is here: https://www.notion.so/cryodrgn-conformational-landscape-analysis-a5af129288d54d1aa95388bdac48235a.
* NEW: Faster training and higher resolution model with Gaussian Fourier featurization (Use `--pe-type gaussian`)
* NEW: `cryodrgn_utils <command> -h` for standalone utility scripts
* NEW: `cryodrgn_utils write_star` for converting cryoDRGN particle selections to `.star` files
* Add pytorch native mixed precision training and fix support for pytorch 1.9+

</details>

<details><summary>Version 0.3.4</summary>

* FIX: Bug in `write_starfile.py` when provided particle stack is chunked (.txt file)
* Support micrograph coordinates and additional column headers to `write_starfile.py`
* New helper scripts: `analyze_convergence.py` (_in beta testing_) contributed by <a href="bmp@mit.edu">Barrett Powell</a> (thanks!) and `make_random_selection.py` for splitting up particle stacks for training

</details>

<details><summary>Version 0.3.3</summary>

* Faster image preprocessing and smaller memory footprint
* New: `cryodrgn preprocess` for large datasets (_in beta testing_ - see <a href="https://www.notion.so/cryodrgn-preprocess-d84a9d9df8634a6a8bfd32d6b5e737ef">this Notion doc</a> for details)
* Known <a href="https://github.com/zhonge/cryodrgn/issues/66">issue</a> with PyTorch version 1.9+

</details>

<details><summary>Version 0.3.2</summary>

* New: cryoDRGN_filtering.ipynb for interactive filtering/selection of images from the dataset
* New: `cryodrgn view_config`
* Minor performance improvements and compatibility fixes

</details>

<details><summary>Version 0.3.1</summary>

* New: Script `write_starfile.py` to convert (filtered) particle selection to a .star file
* More visualizations in `cryodrgn analyze`

</details>

<details><summary>Version 0.3.0</summary>

* New: GPU parallelization with flag `--multigpu`
* New: Mode for accelerated mixed precision training with flag `--amp`, available for NVIDIA tensor core GPUs
* Interface update:
    * Renamed encoder arguments `--qdim` and `--qlayers` to `--enc-dim` and `--enc-layers`
    * Renamed decoder arguments `--pdim` and `--players` to `--dec-dim` and `--dec-layers`
* Argument default changes:
    * Flipped the default for `--invert-data` to True by default
    * Flipped the default for `--window` to True by default
* Updated training recommendations in below quick start guide
* Updates to cryodrgn analyze
    * More visualizations
    * Order kmeans volumes according to distances in latent space (previously random)
    * More features for particle selection and filtering in the Jupiter notebook

</details>


<details><summary>Version 0.2.1</summary>

* New: Parsing of RELION 3.1 files
* Fix: Compatibility with pytorch 1.5

</details>


<details><summary>Version 0.2.0</summary>

* New interface and proper python packaging with `setup.py`. This version has identical functionality and argument usage as previous versions, however tools are now available from a common entry point. See:

    `$ cryodrgn <command> -h`

* New analysis pipeline `cryodrgn analyze`
* New latent space traversal scripts with `cryodrgn graph_traversal` and `cryodrgn pc_traversal`.

</details>

## Installation:

`cryodrgn` may be installed via `pip`, and we recommend installing `cryodrgn` in a clean conda environment.

    # Create and activate conda environment
    (base) $ conda create --name cryodrgn python=3.9
    (cryodrgn) $ conda activate cryodrgn

    # install cryodrgn
    (cryodrgn) $ pip install cryodrgn

More installation instructions are found in the [documentation](https://ez-lab.gitbook.io/cryodrgn/installation).

## Quickstart: heterogeneous reconstruction with consensus poses

### 1. Preprocess image stack

First resize your particle images using the `cryodrgn downsample` command:

<details><summary><code>$ cryodrgn downsample -h</code></summary>

    usage: cryodrgn downsample [-h] -D D -o MRCS [--is-vol] [--chunk CHUNK]
                               [--datadir DATADIR]
                               mrcs

    Downsample an image stack or volume by clipping fourier frequencies

    positional arguments:
      mrcs               Input images or volume (.mrc, .mrcs, .star, .cs, or .txt)

    optional arguments:
      -h, --help         show this help message and exit
      -D D               New box size in pixels, must be even
      -o MRCS            Output image stack (.mrcs) or volume (.mrc)
      --is-vol           Flag if input .mrc is a volume
      --chunk CHUNK      Chunksize (in # of images) to split particle stack when
                         saving
      --relion31         Flag for relion3.1 star format
      --datadir DATADIR  Optionally provide path to input .mrcs if loading from a
                         .star or .cs file
      --max-threads MAX_THREADS
                         Maximum number of CPU cores for parallelization (default: 16)
      --ind PKL          Filter image stack by these indices

</details>

We recommend first downsampling images to 128x128 since larger images can take much longer to train:

    $ cryodrgn downsample [input particle stack] -D 128 -o particles.128.mrcs

The maximum recommended image size is D=256, so we also recommend downsampling your images to D=256 if your images are larger than 256x256:

    $ cryodrgn downsample [input particle stack] -D 256 -o particles.256.mrcs

The input file format can be a single `.mrcs` file, a `.txt` file containing paths to multiple `.mrcs` files, a RELION `.star` file, or a cryoSPARC `.cs` file. For the latter two options, if the relative paths to the `.mrcs` are broken, the argument `--datadir` can be used to supply the path to where the `.mrcs` files are located.

If there are memory issues with downsampling large particle stacks, add the `--chunk 10000` argument to save images as separate `.mrcs` files of 10k images.

### 2. Parse image poses from a consensus homogeneous reconstruction

CryoDRGN expects image poses to be stored in a binary pickle format (`.pkl`). Use the `parse_pose_star` or `parse_pose_csparc` command to extract the poses from a `.star` file or a `.cs` file, respectively.

Example usage to parse image poses from a RELION 3.1 starfile:

    $ cryodrgn parse_pose_star particles.star -o pose.pkl -D 300

Example usage to parse image poses from a cryoSPARC homogeneous refinement particles.cs file:

    $ cryodrgn parse_pose_csparc cryosparc_P27_J3_005_particles.cs -o pose.pkl -D 300

**Note:** The `-D` argument should be the box size of the consensus refinement (and not the downsampled images from step 1) so that the units for translation shifts are parsed correctly.

### 3. Parse CTF parameters from a .star/.cs file

CryoDRGN expects CTF parameters to be stored in a binary pickle format (`.pkl`). Use the `parse_ctf_star` or `parse_ctf_csparc` command to extract the relevant CTF parameters from a `.star` file or a `.cs` file, respectively.

Example usage for a .star file:

    $ cryodrgn parse_ctf_star particles.star -D 300 --Apix 1.03 -o ctf.pkl

The `-D` and `--Apix` arguments should be set to the box size and Angstrom/pixel of the original `.mrcs` file (before any downsampling).

Example usage for a .cs file:

    $ cryodrgn parse_ctf_csparc cryosparc_P27_J3_005_particles.cs -o ctf.pkl

### 4. (Optional) Test pose/CTF parameters parsing

Next, test that pose and CTF parameters were parsed correctly using the voxel-based backprojection script.
The goal is to quickly verify that there are no major problems with the extracted values and that the output structure resembles the structure from the consensus reconstruction before training.

Example usage:

    $ cryodrgn backproject_voxel projections.128.mrcs \
            --poses pose.pkl \
            --ctf ctf.pkl \
            -o backproject.128.mrc

The output structure `backproject.128.mrc` will not match the consensus reconstruction exactly as the `backproject_voxel` command backprojects phase-flipped particles onto the voxel grid, and by default only uses the first 10k images. If the structure is too noisy, you can increase the number of images that are used with the `--first` argument.

**Note:** If the volume does not resemble your structure, you may need to use the flag `--uninvert-data`. This flips the data sign (e.g. light-on-dark or dark-on-light), which may be needed depending on the convention used in upstream processing tools.

### 5. Running cryoDRGN heterogeneous reconstruction

When the input images (.mrcs), poses (.pkl), and CTF parameters (.pkl) have been prepared, a cryoDRGN model can be trained with following command:

<details><summary><code>$ cryodrgn train_vae -h</code></summary>

	usage: cryodrgn train_vae [-h] -o OUTDIR --zdim ZDIM --poses POSES [--ctf pkl]
	                          [--load WEIGHTS.PKL] [--checkpoint CHECKPOINT]
	                          [--log-interval LOG_INTERVAL] [-v] [--seed SEED]
	                          [--ind PKL] [--uninvert-data] [--no-window]
	                          [--window-r WINDOW_R] [--datadir DATADIR] [--lazy]
	                          [--preprocessed] [--max-threads MAX_THREADS]
	                          [--tilt TILT] [--tilt-deg TILT_DEG] [-n NUM_EPOCHS]
	                          [-b BATCH_SIZE] [--wd WD] [--lr LR] [--beta BETA]
	                          [--beta-control BETA_CONTROL] [--norm NORM NORM]
	                          [--no-amp] [--multigpu] [--do-pose-sgd]
	                          [--pretrain PRETRAIN] [--emb-type {s2s2,quat}]
	                          [--pose-lr POSE_LR] [--enc-layers QLAYERS]
	                          [--enc-dim QDIM]
	                          [--encode-mode {conv,resid,mlp,tilt}]
	                          [--enc-mask ENC_MASK] [--use-real]
	                          [--dec-layers PLAYERS] [--dec-dim PDIM]
	                          [--pe-type {geom_ft,geom_full,geom_lowf,geom_nohighf,linear_lowf,gaussian,none}]
	                          [--feat-sigma FEAT_SIGMA] [--pe-dim PE_DIM]
	                          [--domain {hartley,fourier}]
	                          [--activation {relu,leaky_relu}]
	                          particles

	Train a VAE for heterogeneous reconstruction with known pose

	positional arguments:
	  particles             Input particles (.mrcs, .star, .cs, or .txt)

	optional arguments:
	  -h, --help            show this help message and exit
	  -o OUTDIR, --outdir OUTDIR
	                        Output directory to save model
	  --zdim ZDIM           Dimension of latent variable
	  --poses POSES         Image poses (.pkl)
	  --ctf pkl             CTF parameters (.pkl)
	  --load WEIGHTS.PKL    Initialize training from a checkpoint
	  --checkpoint CHECKPOINT
	                        Checkpointing interval in N_EPOCHS (default: 1)
	  --log-interval LOG_INTERVAL
	                        Logging interval in N_IMGS (default: 1000)
	  -v, --verbose         Increaes verbosity
	  --seed SEED           Random seed

	Dataset loading:
	  --ind PKL             Filter particle stack by these indices
	  --uninvert-data       Do not invert data sign
	  --no-window           Turn off real space windowing of dataset
	  --window-r WINDOW_R   Windowing radius (default: 0.85)
	  --datadir DATADIR     Path prefix to particle stack if loading relative
	                        paths from a .star or .cs file
	  --lazy                Lazy loading if full dataset is too large to fit in
	                        memory (Should copy dataset to SSD)
	  --preprocessed        Skip preprocessing steps if input data is from
	                        cryodrgn preprocess_mrcs
	  --max-threads MAX_THREADS
	                        Maximum number of CPU cores for FFT parallelization
	                        (default: 16)

	Tilt series:
	  --tilt TILT           Particles (.mrcs)
	  --tilt-deg TILT_DEG   X-axis tilt offset in degrees (default: 45)

	Training parameters:
	  -n NUM_EPOCHS, --num-epochs NUM_EPOCHS
	                        Number of training epochs (default: 20)
	  -b BATCH_SIZE, --batch-size BATCH_SIZE
	                        Minibatch size (default: 8)
	  --wd WD               Weight decay in Adam optimizer (default: 0)
	  --lr LR               Learning rate in Adam optimizer (default: 0.0001)
	  --beta BETA           Choice of beta schedule or a constant for KLD weight
	                        (default: 1/zdim)
	  --beta-control BETA_CONTROL
	                        KL-Controlled VAE gamma. Beta is KL target. (default:
	                        None)
	  --norm NORM NORM      Data normalization as shift, 1/scale (default: 0, std
	                        of dataset)
	  --no-amp              Do not use mixed-precision training
	  --multigpu            Parallelize training across all detected GPUs

	Pose SGD:
	  --do-pose-sgd         Refine poses with gradient descent
	  --pretrain PRETRAIN   Number of epochs with fixed poses before pose SGD
	                        (default: 1)
	  --emb-type {s2s2,quat}
	                        SO(3) embedding type for pose SGD (default: quat)
	  --pose-lr POSE_LR     Learning rate for pose optimizer (default: 0.0003)

	Encoder Network:
	  --enc-layers QLAYERS  Number of hidden layers (default: 3)
	  --enc-dim QDIM        Number of nodes in hidden layers (default: 1024)
	  --encode-mode {conv,resid,mlp,tilt}
	                        Type of encoder network (default: resid)
	  --enc-mask ENC_MASK   Circular mask of image for encoder (default: D/2; -1
	                        for no mask)
	  --use-real            Use real space image for encoder (for convolutional
	                        encoder)

	Decoder Network:
	  --dec-layers PLAYERS  Number of hidden layers (default: 3)
	  --dec-dim PDIM        Number of nodes in hidden layers (default: 1024)
	  --pe-type {geom_ft,geom_full,geom_lowf,geom_nohighf,linear_lowf,gaussian,none}
	                        Type of positional encoding (default: gaussian)
	  --feat-sigma FEAT_SIGMA
	                        Scale for random Gaussian features
	  --pe-dim PE_DIM       Num features in positional encoding (default: image D)
	  --domain {hartley,fourier}
	                        Decoder representation domain (default: fourier)
	  --activation {relu,leaky_relu}
	                        Activation (default: relu)

</details>

Many of the parameters of this script have sensible defaults. The required arguments are:

* an input image stack (`.mrcs` or other listed file types)
* `--poses`, image poses (`.pkl`) that correspond to the input images
* `--ctf`, ctf parameters (`.pkl`), unless phase-flipped images are used
* `--zdim`, the dimension of the latent variable
* `-o`, a clean output directory for saving results

Additional parameters which are typically set include:

* `-n`, Number of epochs to train
* `--uninvert-data`, Use if particles are dark on light (negative stain format)
* Architecture parameters with `--enc-layers`, `--enc-dim`, `--dec-layers`, `--dec-dim`
* `--multigpu` to enable parallelized training across multiple GPUs

### Recommended usage:

1) It is highly recommended to first train on lower resolution images (e.g. D=128) to sanity check results and perform any particle filtering.

Example command to train a cryoDRGN model for 25 epochs on an image dataset `projections.128.mrcs` with poses `pose.pkl` and ctf parameters `ctf.pkl`:

    # 8-D latent variable model, small images
    $ cryodrgn train_vae projections.128.mrcs \
            --poses pose.pkl \
            --ctf ctf.pkl \
            --zdim 8 -n 25 \
            -o 00_cryodrgn128

2) After validation, pose optimization, and any necessary particle filtering, then train on the full resolution images (up to D=256):

Example command to train a cryoDRGN model for 25 epochs on an image dataset `projections.256.mrcs` with poses `pose.pkl` and ctf parameters `ctf.pkl`:

    # 8-D latent variable model, larger images
    $ cryodrgn train_vae projections.256.mrcs \
            --poses pose.pkl \
            --ctf ctf.pkl \
            --zdim 8 -n 25 \
            -o 01_cryodrgn256

The number of epochs `-n` refers to the number of full passes through the dataset for training, and should be modified depending on the number of particles in the dataset. For a 100k particle dataset on 1 V100 GPU, the above settings required ~12 min/epoch for D=128 images and ~47 min/epoch for D=256 images.

If you would like to train longer, a training job can be extended with the `--load` argument. For example to extend the training of the previous example to 50 epochs:

    $ cryodrgn train_vae projections.256.mrcs \
            --poses pose.pkl \
            --ctf ctf.pkl \
            --zdim 8 -n 50 \
            -o 01_cryodrgn256 \
            --load 01_cryodrgn256/weights.24.pkl # 0-based indexing

### Accelerated training with GPU parallelization

Use cryoDRGN's `--multigpu` flag to enable parallelized training across all detected GPUs on the machine. To select specific GPUs for cryoDRGN to run on, use the environmental variable `CUDA_VISIBLE_DEVICES`, e.g.:

    $ cryodrgn train_vae ... # Run on GPU 0
    $ cryodrgn train_vae ... --multigpu # Run on all GPUs on the machine
    $ CUDA_VISIBLE_DEVICES=0,3 cryodrgn train_vae ... --multigpu # Run on GPU 0,3

We recommend using `--multigpu` for large images, e.g. D=256. Note that GPU computation may not be the training bottleneck for smaller images (D=128). In this case, `--multigpu` may not speed up training (while taking up additional compute resources).

With `--multigpu`, the batch size is multiplied by the number of available GPUs to better utilize GPU resources. We note that GPU utilization may be further improved by increasing the batch size (e.g. `-b 16`), however, faster wall-clock time per epoch does not necessarily lead to faster *model training* since the training dynamics are affected (fewer model updates per epoch with larger `-b`), and using `--multigpu` may require increasing the total number of epochs.

### Local pose refinement -- *beta*

Depending on the quality of the consensus reconstruction, image poses may contain errors.
Image poses may be *locally* refined using the `--do-pose-sgd` flag, however, we recommend reaching out to the developers for recommended training settings.

## 6. Analysis of results

Once the model has finished training, the output directory will contain a configuration file `config.yaml`, neural network weights `weights.pkl`, image poses (if performing pose sgd) `pose.pkl`, and the latent embeddings for each image `z.pkl`. The latent embeddings are provided in the same order as the input particles. To analyze these results, use the `cryodrgn analyze` command to visualize the latent space and generate structures. `cryodrgn analyze` will also provide a template jupyter notebook for further interactive visualization and analysis.


### cryodrgn analyze

<details><summary><code>$ cryodrgn analyze -h</code></summary>

	usage: cryodrgn analyze [-h] [--device DEVICE] [-o OUTDIR] [--skip-vol]
	                        [--skip-umap] [--Apix APIX] [--flip] [--invert]
	                        [-d DOWNSAMPLE] [--pc PC] [--ksample KSAMPLE]
	                        workdir epoch

	Visualize latent space and generate volumes

	positional arguments:
	  workdir               Directory with cryoDRGN results
	  epoch                 Epoch number N to analyze (0-based indexing,
	                        corresponding to z.N.pkl, weights.N.pkl)

	optional arguments:
	  -h, --help            show this help message and exit
	  --device DEVICE       Optionally specify CUDA device
	  -o OUTDIR, --outdir OUTDIR
	                        Output directory for analysis results (default:
	                        [workdir]/analyze.[epoch])
	  --skip-vol            Skip generation of volumes
	  --skip-umap           Skip running UMAP

	Extra arguments for volume generation:
	  --Apix APIX           Pixel size to add to .mrc header (default: 1 A/pix)
	  --flip                Flip handedness of output volumes
	  --invert              Invert contrast of output volumes
	  -d DOWNSAMPLE, --downsample DOWNSAMPLE
	                        Downsample volumes to this box size (pixels)
	  --pc PC               Number of principal component traversals to generate
	                        (default: 2)
	  --ksample KSAMPLE     Number of kmeans samples to generate (default: 20)

</details>

This script runs a series of standard analyses:

* PCA visualization of the latent embeddings
* UMAP visualization of the latent embeddings
* Generation of volumes. See note [1].
* Generation of trajectories along the first and second principal components of the latent embeddings
* Generation of template jupyter notebooks that may be used for further interactive analyses, visualization, and volume generation

Example usage to analyze results from the direction `01_cryodrgn256` containing results after 25 epochs of training:

    $ cryodrgn analyze 01_cryodrgn256 24 --Apix 1.31 # 24 for 0-based indexing of epoch numbers

Notes:

[1] Volumes are generated after k-means clustering of the latent embeddings with k=20 by default. Note that we use k-means clustering here not to identify clusters, but to segment the latent space and generate structures from different regions of the latent space. The number of structures that are generated may be increased with the option `--ksample`.

[2] The `cryodrgn analyze` command chains together a series of calls to `cryodrgn eval_vol` and other scripts that can be run separately for more flexibility. These scripts are located in the `analysis_scripts` directory within the source code.

### Generating additional volumes

A simple way of generating additional volumes is to increase the number of k-means samples in `cryodrgn analyze` by using the flag `--ksample 100` (for 100 structures). For additional flexibility, `cryodrgn eval_vol` may be called directly:

<details><summary><code>$ cryodrgn eval_vol -h</code></summary>

	usage: cryodrgn eval_vol [-h] -c PKL -o O [--prefix PREFIX] [-v]
	                         [-z [Z [Z ...]]] [--z-start [Z_START [Z_START ...]]]
	                         [--z-end [Z_END [Z_END ...]]] [-n N] [--zfile ZFILE]
	                         [--Apix APIX] [--flip] [-d DOWNSAMPLE]
	                         [--norm NORM NORM] [-D D] [--enc-layers QLAYERS]
	                         [--enc-dim QDIM] [--zdim ZDIM]
	                         [--encode-mode {conv,resid,mlp,tilt}]
	                         [--dec-layers PLAYERS] [--dec-dim PDIM]
	                         [--enc-mask ENC_MASK]
	                         [--pe-type {geom_ft,geom_full,geom_lowf,geom_nohighf,linear_lowf,none}]
	                         [--pe-dim PE_DIM] [--domain {hartley,fourier}]
	                         [--l-extent L_EXTENT]
	                         [--activation {relu,leaky_relu}]
	                         weights

	Evaluate the decoder at specified values of z

	positional arguments:
	  weights               Model weights

	optional arguments:
	  -h, --help             show this help message and exit
	  -c YAML, --config YAML CryoDRGN config.yaml file
	  -o O                   Output .mrc or directory
	  --prefix PREFIX        Prefix when writing out multiple .mrc files (default: vol_)
	  -v, --verbose          Increase verbosity

	Specify z values:
	  -z [Z [Z ...]]        Specify one z-value
	  --z-start [Z_START [Z_START ...]]
	                        Specify a starting z-value
	  --z-end [Z_END [Z_END ...]]
	                        Specify an ending z-value
	  -n N                  Number of structures between [z_start, z_end]
	  --zfile ZFILE         Text file with z-values to evaluate

	Volume arguments:
	  --Apix APIX           Pixel size to add to .mrc header (default: 1 A/pix)
	  --flip                Flip handedness of output volume
	  -d DOWNSAMPLE, --downsample DOWNSAMPLE
	                        Downsample volumes to this box size (pixels)

	Overwrite architecture hyperparameters in config.yaml:
	  --norm NORM NORM
	  -D D                  Box size
	  --enc-layers QLAYERS  Number of hidden layers
	  --enc-dim QDIM        Number of nodes in hidden layers
	  --zdim ZDIM           Dimension of latent variable
	  --encode-mode {conv,resid,mlp,tilt}
	                        Type of encoder network
	  --dec-layers PLAYERS  Number of hidden layers
	  --dec-dim PDIM        Number of nodes in hidden layers
	  --enc-mask ENC_MASK   Circular mask radius for image encoder
	  --pe-type {geom_ft,geom_full,geom_lowf,geom_nohighf,linear_lowf,none}
	                        Type of positional encoding
	  --pe-dim PE_DIM       Num sinusoid features in positional encoding (default:
	                        D/2)
	  --domain {hartley,fourier}
	  --l-extent L_EXTENT   Coordinate lattice size
	  --activation {relu,leaky_relu}
	                        Activation (default: relu)

</details>

**Example usage:**

To generate a volume at a single value of the latent variable:

    $ cryodrgn eval_vol [YOUR_WORKDIR]/weights.pkl --config [YOUR_WORKDIR]/config.yaml -z ZVALUE -o reconstruct.mrc

The number of inputs for `-z` must match the dimension of your latent variable.

Or to generate a trajectory of structures from a defined start and ending point, use the `--z-start` and `--z-end` arugments:

    $ cryodrgn eval_vol [YOUR_WORKDIR]/weights.pkl --config [YOUR_WORKDIR]/config.yaml --z-start -3 --z-end 3 -n 20 -o [WORKDIR]/trajectory

This example generates 20 structures at evenly spaced values between z=[-3,3], assuming a 1-dimensional latent variable model.

Finally, a series of structures can be generated using values of z given in a file specified by the arugment `--zfile`:

    $ cryodrgn eval_vol [WORKDIR]/weights.pkl --config [WORKDIR]/config.yaml --zfile zvalues.txt -o [WORKDIR]/trajectory

The input to `--zfile` is expected to be an array of dimension (N_volumes x zdim), loaded with np.loadtxt.

### Making trajectories

Two additional commands can be used in conjunction with `cryodrgn eval_vol` to generate trajectories:

    $ cryodrgn pc_traversal -h
    $ cryodrgn graph_traversal -h

These scripts produce a text file of z values that can be input to `cryodrgn eval_vol` to generate a series of structures that can be visualized as a trajectory in ChimeraX (https://www.cgl.ucsf.edu/chimerax).

Documentation: https://ez-lab.gitbook.io/cryodrgn/cryodrgn-graph-traversal-for-making-long-trajectories

### cryodrgn analyze_landscape

NEW in version 1.0: There are two additional tools `cryodrgn analyze_landscape` and `cryodrgn analyze_landscape_full` for more comprehensive and automated analyses of cryodrgn results.

Documentation: https://ez-lab.gitbook.io/cryodrgn/cryodrgn-conformational-landscape-analysis

## CryoDRGN2 for *Ab Initio* Reconstruction

To perform *ab initio* heterogeneous reconstruction, use `cryodrgn abinit_het`. The arguments are similar to `cryodrgn train_vae`, but the `--poses` argument is not required.

For homogeneous reconstruction, use `cryodrgn abinit_homo`.

Documentation: https://ez-lab.gitbook.io/cryodrgn/cryodrgn2-ab-initio-reconstruction

## CryoDRGN-ET for subtomogram analysis

Coming soon in version 3.0.

## References:

For a complete description of the method, see:

* CryoDRGN: reconstruction of heterogeneous cryo-EM structures using neural networks
Ellen D. Zhong, Tristan Bepler, Bonnie Berger*, Joseph H Davis*
Nature Methods 2021, https://doi.org/10.1038/s41592-020-01049-4 [pdf](https://ezlab.princeton.edu/assets/pdf/2021_cryodrgn_nature_methods.pdf)

An earlier version of this work appeared at ICLR 2020:

* Reconstructing continuous distributions of protein structure from cryo-EM images
Ellen D. Zhong, Tristan Bepler, Joseph H. Davis*, Bonnie Berger*
ICLR 2020, Spotlight, https://arxiv.org/abs/1909.05215

CryoDRGN2's ab initio reconstruction algorithms were published at ICCV:

* CryoDRGN2: Ab Initio Neural Reconstruction of 3D Protein Structures From Real Cryo-EM Images
Ellen D. Zhong, Adam Lerer, Joseph H Davis, and Bonnie Berger
International Conference on Computer Vision (ICCV) 2021, [paper](https://openaccess.thecvf.com/content/ICCV2021/papers/Zhong_CryoDRGN2_Ab_Initio_Neural_Reconstruction_of_3D_Protein_Structures_From_ICCV_2021_paper.pdf)

A protocols paper that describes the analysis of the EMPIAR-10076 assembling ribosome dataset:

* Uncovering structural ensembles from single particle cryo-EM data using cryoDRGN
Laurel Kinman, Barrett Powell, Ellen D. Zhong*, Bonnie Berger*, Joseph H Davis*
Nature Protocols 2023, https://doi.org/10.1038/s41596-022-00763-x

## Contact

Please submit any bug reports, feature requests, or general usage feedback as a github issue, or post in the Google Group: https://groups.google.com/g/cryodrgn.
