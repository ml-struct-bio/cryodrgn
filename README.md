# :snowflake::dragon: cryoDRGN: Deep Reconstructing Generative Networks for cryo-EM heterogeneous reconstruction

CryoDRGN is a neural network based algorithm for heterogeneous cryo-EM reconstruction. In particular, the method models a *continuous* distribution over 3D structures by using a neural network based representation for the volume.

## Manuscript:

CryoDRGN: Reconstruction of heterogeneous structures from cryo-electron micrographs using neural networks.
Ellen D. Zhong, Tristan Bepler, Bonnie Berger*, Joseph H. Davis*.
https://www.biorxiv.org/content/10.1101/2020.03.27.003871v1

Reconstructing continuous distributions of 3D protein structure from cryo-EM images.
Ellen D. Zhong, Tristan Bepler, Joseph H. Davis*, Bonnie Berger*.
ICLR 2020, https://arxiv.org/abs/1909.05215

## New in v0.2

**Version 0.2.1:**
* New: Parsing of RELION 3.1 files
* Fix: Compatibility with pytorch 1.5

**Version 0.2.0:**
* New interface and proper python packaing with setup.py. This version has identical functionality and argument usage as previous versions, however tools are now available from a common entry point. See:
    
    `$ cryodrgn <command> -h`

* New analysis pipeline `cryodrgn analyze`
* New latent space traversal scripts with `cryodrgn graph_traversal` and `cryodrgn pc_traversal`.

*Last updated 4/20/2020.*

## Installation/dependencies:

Until the cryoDRGN conda/pip package is available, for now, git clone the source code and install the following dependencies with anaconda, replacing the cudatoolkit version as necessary:

    # Create conda environment
    conda create --name cryodrgn python=3.7
    conda activate cryodrgn

    # Install dependencies
    conda install pytorch cudatoolkit=10.1 -c pytorch # Replace cudatoolkit version if needed
    conda install pandas
    
    # Install dependencies for latent space visualization 
    conda install seaborn scikit-learn 
    conda install -c conda-forge umap-learn
    conda install -c conda-forge jupyterlab
    pip install ipywidgets
    pip install cufflinks

    # Clone source code and install
    git clone https://github.com/zhonge/cryodrgn.git
    cd cryodrgn
    git checkout 0.2.1
    python setup.py install

## Quickstart: heterogeneous reconstruction with consensus poses

### 1. Preprocess image stack

First resize your particle images for initial pilot experiments with cryoDRGN using the `cryodrgn downsample` command:

    $ cryodrgn downsample -h
    usage: cryodrgn downsample [-h] -D D -o MRCS [--is-vol] [--chunk CHUNK]
                               [--datadir DATADIR]
                               mrcs
    
    Downsample an image stack or volume by clipping fourier frequencies
    
    positional arguments:
      mrcs               Input particles or volume (.mrc, .mrcs, .star, or .txt)
    
    optional arguments:
      -h, --help         show this help message and exit
      -D D               New box size in pixels, must be even
      -o MRCS            Output projection stack (.mrcs)
      --is-vol           Flag if input .mrc is a volume
      --chunk CHUNK      Chunksize (in # of images) to split particle stack when
                         saving
      --datadir DATADIR  Optionally provide path to input .mrcs if loading from a
                         .star or .cs file

Since larger images require longer training times, it is recommended to first train cryoDRGN on images downsized to D=128, (128x128 pixel image):
    
    $ cryodrgn downsample [input particle stack] -D 128 -o particles.128.mrcs

The maximum recommended image size is D=256, so it is also recommended to downsample your images to D=256 if your images are larger than 256x256:

    $ cryodrgn downsample [input particle stack] -D 256 -o particles.256.mrcs

The input file format can be a single `.mrcs` file, a `.txt` file containing paths to multiple `.mrcs` files, a `.star` file, or a cryoSPARC `.cs` file. For the latter two options, if the relative paths to the `.mrcs` are broken, the argument `--datadir` can be used to supply the path to where the `.mrcs` files are located. 

If there are memory issues with large particle stacks, add the `--chunk 10000` argument to save out images as separate `.mrcs` files of 10k images. 

### 2. Parse image poses from a consensus homogeneous reconstruction

CryoDRGN expects image poses in a binary pickle format (`.pkl`). Use the `parse_pose_star` or `parse_pose_csparc` command to extract the poses from a `.star` file or a `.cs` file, respectively.

Example usage to parse image poses from a RELION starfile:
    
    $ cryodrgn parse_pose_star particles.star -o pose.pkl -D 300

Example usage to parse image poses from a cryoSPARC homogeneous refinement particles.cs file:

    $ cryodrgn parse_pose_csparc cryosparc_P27_J3_005_particles.cs -o pose.pkl -D 300

The `-D` argument should be set to the box size of the original reconstruction (before any downsampling). 

### 3. Parse CTF parameters from a .star/.cs file

CryoDRGN expects CTF parameters in a binary pickle format (`.pkl`). Use the `parse_ctf_star` or `parse_ctf_csparc` command to extract the relevant CTF parameters from a `.star` file or a `.cs` file, respectively.

Example usage for a .star file:
    
    $ cryodrgn parse_ctf_star particles.star -D 300 --Apix 1.03 -o ctf.pkl

The `-D` and `--Apix` arguments should be set to the box size and Angstrom/pixel of the original `.mrcs` file (before any downsampling). 

Example usage for a .cs file:

    $ cryodrgn parse_ctf_csparc cryosparc_P27_J3_005_particles.cs -o ctf.pkl

### 4. Test pose/CTF parameters parsing

Next, test that pose and CTF parameters were parsed correctly using the voxel-based backprojection script.
The goal is to quickly verify that there are no major problems with the extracted values and that the output structure resembles the structure from the consensus reconstruction before beginning training.

Example usage: 

    $ cryodrgn backproject_voxel projections.128.mrcs \
            --poses pose.pkl \
            --ctf ctf.pkl \
            --invert-data \ # Invert sign of dataset; Use if particles are white on black
            -o backproject.128.mrc

The output structure `backproject.128.mrc` will not match the consensus reconstruction exactly as the `backproject_voxel` command backprojects phase-flipped particles onto the voxel grid, and by default only uses the first 10k images. If the structure is too noisy, you can increase the number of images that are used with the `--first` argument.

### 5. Running cryoDRGN heterogeneous reconstruction

When the input image stack (.mrcs), image poses (.pkl), and CTF parameters (.pkl) have been prepared, a cryoDRGN model can be trained with following script:

    $ cryodrgn train_vae -h
    usage: cryodrgn train_vae [-h] -o OUTDIR --zdim ZDIM --poses POSES [--ctf pkl]
                              [--load LOAD] [--checkpoint CHECKPOINT]
                              [--log-interval LOG_INTERVAL] [-v] [--seed SEED]
                              [--invert-data] [--window] [--ind IND] [--lazy]
                              [--datadir DATADIR] [--tilt TILT]
                              [--tilt-deg TILT_DEG] [-n NUM_EPOCHS]
                              [-b BATCH_SIZE] [--wd WD] [--lr LR] [--beta BETA]
                              [--beta-control BETA_CONTROL] [--norm NORM NORM]
                              [--do-pose-sgd] [--pretrain PRETRAIN]
                              [--emb-type {s2s2,quat}] [--pose-lr POSE_LR]
                              [--qlayers QLAYERS] [--qdim QDIM]
                              [--encode-mode {conv,resid,mlp,tilt}]
                              [--enc-mask ENC_MASK] [--use-real]
                              [--players PLAYERS] [--pdim PDIM]
                              [--pe-type {geom_ft,geom_full,geom_lowf,geom_nohighf,linear_lowf,none}]
                              [--domain {hartley,fourier}]
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
      --ctf pkl             CTF parameters (.pkl) if particle stack is not phase
                            flipped
      --load LOAD           Initialize training from a checkpoint
      --checkpoint CHECKPOINT
                            Checkpointing interval in N_EPOCHS (default: 1)
      --log-interval LOG_INTERVAL
                            Logging interval in N_IMGS (default: 1000)
      -v, --verbose         Increaes verbosity
      --seed SEED           Random seed
    
    Dataset loading:
      --invert-data         Invert data sign
      --window              Real space windowing of dataset
      --ind IND             Filter particle stack by these indices
      --lazy                Lazy loading if full dataset is too large to fit in
                            memory
      --datadir DATADIR     Path prefix to particle stack if loading relative
                            paths from a .star or .cs file
    
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
                            (default: 1.0)
      --beta-control BETA_CONTROL
                            KL-Controlled VAE gamma. Beta is KL target. (default:
                            None)
      --norm NORM NORM      Data normalization as shift, 1/scale (default: mean,
                            std of dataset)
      --do-pose-sgd         Refine poses with gradient descent
      --pretrain PRETRAIN   Number of epochs with fixed poses before pose SGD
                            (default: 1)
      --emb-type {s2s2,quat}
                            SO(3) embedding type for pose SGD (default: quat)
      --pose-lr POSE_LR     Learning rate for pose optimizer (default: 0.0003)
    
    Encoder Network:
      --qlayers QLAYERS     Number of hidden layers (default: 3)
      --qdim QDIM           Number of nodes in hidden layers (default: 256)
      --encode-mode {conv,resid,mlp,tilt}
                            Type of encoder network (default: resid)
      --enc-mask ENC_MASK   Circular mask of image for encoder (default: D/2; -1
                            for no mask)
      --use-real            Use real space image for encoder (for convolutional
                            encoder)
    
    Decoder Network:
      --players PLAYERS     Number of hidden layers (default: 3)
      --pdim PDIM           Number of nodes in hidden layers (default: 256)
      --pe-type {geom_ft,geom_full,geom_lowf,geom_nohighf,linear_lowf,none}
                            Type of positional encoding (default: geom_lowf)
      --domain {hartley,fourier}
                            Decoder representation domain (default: fourier)
    
Many of the parameters of this script have sensisible defaults. The required arguments are:

* an input image stack (`.mrcs` or other listed file types)
* `--poses`, image poses (`.pkl`) which correspond to the input images
* `--zdim`, the dimension of the latent variable
* `-o`, a clean output directory for storing results

Additional parameters which are typically set include:

* `--ctf`, ctf parameters (`.pkl`), unless phase-flipped images are used
* `-n`, Number of epochs to train
* `--invert-data`, Use if particles are white on black
* Architecture parameters with `--qlayers`, `--qdim`, `--players`, `--pdim`

### Example usage:

It is recommended to first train on lower resolution images (e.g. D=128) with `--zdim 1` and with `--zdim 10` using the default architecture (fast). After validation, pose optimization, and any necessary particle filtering, then train on the full resolution image stack (up to D=256) with a large architecture (slow). 

Example command to train a cryoDRGN model for 50 epochs on an image dataset `projections.128.mrcs` with poses `pose.pkl` and ctf parameters `ctf.pkl`:

    # 1-D latent variable model
    $ cryodrgn train_vae projections.128.mrcs 
            --poses pose.pkl \
            --ctf ctf.pkl \
            --zdim 1 -n 50 \
            -o 00_vae128_z1

    # 10-D latent variable model
    $ cryodrgn train_vae projections.128.mrcs 
            --poses pose.pkl \
            --ctf ctf.pkl \
            --zdim 10 -n 50 \
            -o 01_vae128_z10

Example command to train a 10-D latent variable cryoDRGN model for 20 epochs on an image dataset `projections.256.mrcs` with poses `pose.pkl` and ctf parameters `ctf.pkl`:

    # 10-D latent variable model with a larger MLP architecture
    $ cryodrgn train_vae projections.256.mrcs 
            --poses pose.pkl \
            --ctf ctf.pkl \
            --zdim 10 -n 20 \
            --qdim 1024 --qlayers 3 --pdim 1024 --players 3 \
            -o 02_vae256_z10

The number of epochs `-n` refers to the number of full passes through the dataset for training, and should be modified depending on the number of particles in the dataset. For a 100k particle dataset, the above settings required ~6 min per epoch for D=128 images + default architecture, ~12 min/epoch for D=128 images + large architecture, and ~47 min per epoch for D=256 images + large architecture. 

If you would like to train longer, a training job can be extended with the `--load` argument. For example to extend the training of the previous example to 50 epochs:

    $ cryodrgn train_vae projections.256.mrcs
            --poses pose.pkl \
            --ctf ctf.pkl \
            --zdim 10 -n 50 \
            --qdim 1024 --qlayers 3 --pdim 1024 --players 3 \
            -o 02_vae256_z10
            --load 02_vae_256_z10/weights.19.pkl # 0-based indexing

Note: While these settings worked well for the datasets we've tested, they are highly experimental for the general case as different datasets have diverse sources of heterogeneity. Please reach out to the authors with questions/consult -- we'd love to learn more.

### Local pose refinement -- BETA!

Depending on the quality of the consensus reconstruction, image poses may contain errors.
Image poses may be *locally* refined using the `--do-pose-sgd` flag. More details on this method to come!

## 6. Analysis of results

Once the model has finished training, the output directory will contain a configuration file `config.pkl`, neural network weights `weights.pkl`, image poses (if performing pose sgd) `pose.pkl`, and the predicted latent encoding for each image `z.pkl`. To analyze these results, use the `cryodrgn analyze` command to visualize the latent space and generate structures. `cryodrgn analyze` will also provide a template jupyter notebook for further interactive visualization and analysis.

### cryodrgn analyze

    $ cryodrgn analyze -h
    usage: cryodrgn analyze [-h] [--device DEVICE] [--Apix APIX] [--flip]
                            [-d DOWNSAMPLE]
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
    
    Extra arguments for volume generation:
      --Apix APIX           Pixel size to add to .mrc header (default: 1 A/pix)
      --flip                Flip handedness of output volume
      -d DOWNSAMPLE, --downsample DOWNSAMPLE
                            Downsample volumes to this box size (pixels)

This script runs a series of standard analyses:

* PCA of the latent space
* UMAP embedding of the latent space
* Generation of volumes from the latent space. See note [1].
* Generation of trajectories along the first and second principal components
* Generation of a template jupyter notebook that may be used for further interactive analyses and visualization

Example usage to analyze results from the direction `02_vae_256_z10` containing results after 50 epochs of training:

    $ cryodrgn analyze 02_vae_256_z10 49 --Apix 1.31

Notes:

[1] By default, volumes are generated at k-means cluster centers with k=20. Note that we use k-means clustering here not to identify cluters, but to efficiently provide a diverse sample of structures that cover the latent space. For clustering of the latent space, we recommend performing this analysis in the provided jupyter notebook using your favorite clustering algorithm (https://scikit-learn.org/stable/modules/clustering.html).

[2] The `cryodrgn analyze` command chains together a series of calls to `cryodrgn eval_vol` and scripts that can be run separately for more flexibility. The scripts are located in the `analysis_scripts` directory within the source code.

### Generating additional volumes

Additional structures may be generated using `cryodrgn eval_vol`:

    $ cryodrgn eval_vol -h
    usage: cryodrgn eval_vol [-h] [-c PKL] -o O [--prefix PREFIX] [-v]
                             [-z [Z [Z ...]]] [--z-start [Z_START [Z_START ...]]]
                             [--z-end [Z_END [Z_END ...]]] [-n N] [--zfile ZFILE]
                             [--Apix APIX] [--flip] [-d DOWNSAMPLE]
                             [--norm NORM NORM] [-D D] [--qlayers QLAYERS]
                             [--qdim QDIM] [--zdim ZDIM]
                             [--encode-mode {conv,resid,mlp,tilt}]
                             [--players PLAYERS] [--pdim PDIM]
                             [--enc-mask ENC_MASK]
                             [--pe-type {geom_ft,geom_full,geom_lowf,geom_nohighf,linear_lowf,none}]
                             [--domain {hartley,fourier}] [--l-extent L_EXTENT]
                             weights
    
    Evaluate the decoder at specified values of z
    
    positional arguments:
      weights               Model weights
    
    optional arguments:
      -h, --help            show this help message and exit
      -c PKL, --config PKL  CryoDRGN configuration
      -o O                  Output .mrc or directory
      --prefix PREFIX       Prefix when writing out multiple .mrc files (default:
                            vol_)
      -v, --verbose         Increaes verbosity
    
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
    
    Overwrite architecture hyperparameters in config.pkl:
      --norm NORM NORM
      -D D                  Box size
      --qlayers QLAYERS     Number of hidden layers
      --qdim QDIM           Number of nodes in hidden layers
      --zdim ZDIM           Dimension of latent variable
      --encode-mode {conv,resid,mlp,tilt}
                            Type of encoder network
      --players PLAYERS     Number of hidden layers
      --pdim PDIM           Number of nodes in hidden layers
      --enc-mask ENC_MASK   Circular mask radius for image encoder
      --pe-type {geom_ft,geom_full,geom_lowf,geom_nohighf,linear_lowf,none}
                            Type of positional encoding
      --domain {hartley,fourier}
      --l-extent L_EXTENT   Coordinate lattice size

**Example usage:**

To generate a volume at a single value of the latent variable:

    $ cryodrgn eval_vol [YOUR_WORKDIR]/weights.pkl --config [YOUR_WORKDIR]/config.pkl -z ZVALUE -o reconstruct.mrc

The number of inputs for `-z` must match the dimension of your latent variable.

Or to generate a trajectory of structures from a defined start and ending point, use the `--z-start` and `--z-end` arugments:

    $ cryodrgn eval_vol [YOUR_WORKDIR]/weights.pkl --config [YOUR_WORKDIR]/config.pkl --z-start -3 --z-end 3 -n 20 -o [WORKDIR]/trajectory

This example generates 20 structures at evenly spaced values between z=[-3,3]. 

Finally, a series of structures can be generated using values of z given in a file specified by the arugment `--zfile`:

    $ cryodrgn eval_vol [WORKDIR]/weights.pkl --config [WORKDIR]/config.pkl --zfile zvalues.txt -o [WORKDIR]/trajectory

The input to `--zfile` is expected to be an array of dimension (N_volumes x zdim), loaded with np.loadtxt.

### Making trajectories

Two additional commands can be used in conjunction with `cryodrgn eval_vol` to generate trajectories:

    $ cryodrgn pc_traversal -h
    $ cryodrgn graph_traversal -h
    
These scripts provide a text file of z values that can be input to `cryodrgn eval_vol` to generate a series of structures that can be visualized as a trajectory in ChimeraX (https://www.cgl.ucsf.edu/chimerax).

## Fully unsupervised reconstruction

Please reach out to Ellen Zhong (zhonge[at]mit[dot]edu) if you'd like to collaborate on reconstruction of highly heterogeneous datasets where a consensus reconstruction is unavailable.  

## Contact

More documentation and tutorials to come! Bugs reports, feature requests, or general usage feedback to zhonge[at]mit[dot]edu.

