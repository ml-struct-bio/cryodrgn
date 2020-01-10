# :snowflake::dragon: cryoDRGN: Deep Reconstructing Generative Networks for cryo-EM heterogeneous reconstruction

CryoDRGN is a neural network based algorithm for heterogeneous cryo-EM reconstruction. In particular, the method models a *continuous* distribution over 3D structures by using a neural network based representation for the volume.

## Preprint:

Reconstructing continuously heterogeneous structures from single particle cryo-EM with deep generative models.
Ellen D. Zhong, Tristan Bepler, Joseph H. Davis*, Bonnie Berger*
https://arxiv.org/abs/1909.05215

## Installation/dependencies:

Until the cryoDRGN conda package is available, for now, git clone the source code and install the following dependencies with anaconda, replacing the cudatoolkit version as necessary:

    conda create --name cryodrgn
    conda activate cryodrgn
    conda install pytorch=1.0.1 torchvision cudatoolkit=10.0 -c pytorch

Additional requirements for latent space analysis and interactive visualization:

    export CDRGN_SRC="path/to/git/repo"
    conda install seaborn scikit-learn 
    conda install -c conda-forge umap-learn
    conda install -c conda-forge jupyterlab
    pip install ipywidgets
    jupyter nbextension enable --py widgetsnbextension
    pip install cufflinks

## Quickstart: heterogeneous reconstruction with consensus alignments

### 1. Preprocess image stack

Training cryoDRGN networks has not been tested on image sizes above D=256. If your images are larger than D=256, use the following utility to downsample the images:

    $ python $CDRGN_SRC/utils/fouriershrink.py [input particle stack] -D 256 -o [output particle stack] --out-png projections.256.png

It is also recommended to create image stacks at lower resolution (e.g. D=128) for initial testing and pilot experiments with cryoDRGN.

    $ python $CDRGN_SRC/utils/fouriershrink.py [input particle stack] -D 128 -o [output particle stack] --out-png projections.128.png

### 2. Parse alignments from a consensus homogeneous reconstruction

To parse alignments from a RELION starfile:
    
    $ python $CDRGN_SRC/utils/parse_star_alignments.py particles.star -o consensus

To parse alignments from a cryoSPARC homogeneous refinement particles.cs file:

    $ python $CDRGN_SRC/utils/parse_cs_alignments.py cryosparc_P27_J3_005_particles.cs -o consensus

### 3. Parse CTF parameters from a .star file

CryoDRGN currently takes CTF parameters in a binary pickle format (.pkl). Use the utility script `parse_ctf_star.py` to extract the relevant CTF parameters from a .star file. 

Example usage:

    $ python $CDRGN_SRC/utils/parse_ctf_star.py particles.star -N 101845 -D 256 --Apix 1.7 -o ctf.256.pkl

Note: the pixel size is saved in the CTF pickle. Run this script multiple times for each pixel size if cryoDRGN will be run on variable image size particles.

    $ python $CDRGN_SRC/utils/parse_ctf_star.py particles.star -N 101845 -D 128 --Apix 3.4 -o ctf.128.pkl

### 4. Test alignments/CTF parameters were parsed correctly

Test that alignments and CTF parameters were parsed correctly using the voxel-based backprojection script:

    $ python $CDRGN_SRC/backproject_voxel.py -h
    usage: backproject_voxel.py [-h] --poses [POSES [POSES ...]] [--tscale TSCALE]
                                [--ctf pkl] -o O [--invert-data]
                                [--datadir DATADIR] [--ind IND] [--first FIRST]
                                [--tilt TILT] [--tilt-deg TILT_DEG]
                                mrcs
    
    Backproject a stack of images via linear interpolation
    
    positional arguments:
      mrcs                  Input .mrcs image stack
    
    optional arguments:
      -h, --help            show this help message and exit
      --poses [POSES [POSES ...]]
                            Image rotations and optionally translations (.pkl)
      --tscale TSCALE       Scale all translations by this amount (default: 1)
      --ctf pkl             CTF parameters (.pkl) if particle stack is not phase
                            flipped
      -o O                  Output .mrc file
    
    Dataset loading options:
      --invert-data         Invert data sign
      --datadir DATADIR     Path prefix to particle stack if loading relative
                            paths from a .star or .cs file
      --ind IND             Indices to iterate over (pkl)
      --first FIRST         Backproject the first N images (default: 5000)
    
    Tilt series options:
      --tilt TILT           Tilt series .mrcs image stack
      --tilt-deg TILT_DEG   Right-handed x-axis tilt offset in degrees (default:
                            45)
    
Example usage:

    $ python $CDRGN_SRC/backproject_voxel.py projections.128.mrcs \
            --poses consensus.rot.pkl consensus.trans.pkl \
            --tscale 0.4 \
            --ctf ctf.128.pkl \
            --invert-data \ # Invert sign of dataset; Use if particles are white on black
            --first 10000 \
            -o backproject.128.mrc

Check that the output structure `backproject.128.mrc` resembles the structure from the consensus reconstruction. 
It will not match exactly as the `backproject_voxel.py` script performs linear interpolation of phase-flipped particles onto the voxel grid.

NOTE: Parsed translations are given in units of pixels. Use the `--tscale` flag to renormalize the translations if the images you are backprojecting are a different size than the images used in the consensus reconstruction.

### 5. Running cryoDRGN heterogeneous reconstruction

When the input image stack (.mrcs), image poses (.pkl), and CTF parameters (.pkl) have been prepared, a cryoDRGN model can be trained with following script:

    $ python $CDRGN_SRC/vae_het.py -h
    usage: vae_het.py [-h] -o OUTDIR --zdim ZDIM --poses [POSES [POSES ...]]
                      [--tscale TSCALE] [--ctf pkl] [--load LOAD]
                      [--checkpoint CHECKPOINT] [--log-interval LOG_INTERVAL] [-v]
                      [--seed SEED] [--invert-data] [--window] [--ind IND]
                      [--lazy] [--datadir DATADIR] [--tilt TILT]
                      [--tilt-deg TILT_DEG] [-n NUM_EPOCHS] [-b BATCH_SIZE]
                      [--wd WD] [--lr LR] [--beta BETA]
                      [--beta-control BETA_CONTROL] [--norm NORM NORM]
                      [--do-pose-sgd] [--pretrain PRETRAIN]
                      [--emb-type {s2s2,quat}] [--pose-lr POSE_LR]
                      [--qlayers QLAYERS] [--qdim QDIM]
                      [--encode-mode {conv,resid,mlp,tilt}] [--enc-mask ENC_MASK]
                      [--use-real] [--players PLAYERS] [--pdim PDIM]
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
      --poses [POSES [POSES ...]]
                            Image rotations and translations (.pkl)
      --tscale TSCALE       Scale translations by this amount
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
                            Number of training epochs (default: 10)
      -b BATCH_SIZE, --batch-size BATCH_SIZE
                            Minibatch size (default: 10)
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

* an input image stack (.mrcs)
* `--poses`, image poses (.pkl)
* `--zdim`, the dimension of the latent variable
* `-o`, a clean output directory for storing results

Additional parameters which are typically set include:

* `--ctf`, CTF parameters (.pkl), unless phase flipped images are used
* `-n`, Number of epochs to train
* `--invert-data`, depending on the data sign convention used in previous processing steps
* `--t-scale`, for modifying the scale of the translations

### Example usage:

Example command to train a 10-D latent variable cryoDRGN model for 50 epochs on an image dataset `projections.128.mrcs` with poses `consensus.rot.pkl, consensus.trans.pkl` and ctf parameters `ctf.128.pkl`:

    $ python $CDRGN_SRC/vae_het.py projections.128.mrcs \
            --poses consensus.rot.pkl consensus.trans.pkl \
            --tscale .4 \ 
            --ctf ctf.128.pkl \
            --zdim 10 \
            -n 50 \
            -o 00_vae128_z10

* Results will be saved in the specified directory `00_vae128_z10`.
* NOTE: Since translations are stored in units of pixels, `--tscale` may be used to renormalize translations if the consensus poses were obtained with a different box size. 
E.g. Use `--tscale .4` if the consensus reconstruction was performed at a box size of 320 pixels and the cryoDRGN model is trained on 128x128 images.

### Recommended settings:

Pilot experiments:
* Train a model using the default neural network architecture, small image sizes (e.g. D=128 or smaller), and zD=1
* Train a model using the default neural network architecture, small image sizes (e.g. D=128 or smaller), and zD=10

After validation, pose optimization, and any necessary particle filtering, train the full resolution image stack (up to D=256) with a large model:
* Large architecture (e.g. 1024 dim x 3 layer network), large image sizes, zD=10

Note these settings are highly experimental.

### Local pose refinement

Depending on the quality of the consensus reconstruction, image poses may contain errors.
Image poses may be *locally* refined using the `--do-pose-sgd` flag. 

### Analysis

Once the model has finished training, the working directory will contain a configuration file `config.pkl`, neural network weights `weights.pkl`, image poses (if performing pose sgd) `pose.pkl`, and the predicted latent encoding for each image `z.pkl`. The directory will also contain a volume `reconstruct.mrc` evaluated at the mean z value of the dataset. 

To analyze and visualize the learned latent space, use the scripts in the `utils/analysis` subdirectory. Additional structures may be generated using the `eval_decoder.py` script. For example:

    $ python $CDRGN_SRC/eval_decoder.py [WORKDIR]/weights.pkl --config [WORKDIR]/config.pkl -z ZVALUE -o reconstruct.sample1.mrc

Or to generate a trajectory using values of z given in a file `zvalues.txt`:

    $ python $CDRGN_SRC/eval_decoder.py [WORKDIR]/weights.pkl --config [WORKDIR]/config.pkl --zfile zvalues.txt -o [WORKDIR]/trajectory

For models with higher dimensional latent variables (zD>2), a shell script is provided which samples 20 structures from the latent space, runs UMAP dimensionality reduction, and creates a template jupyter-notebook in the working directory, which may be used for interactive visualization of the results:

    $ $CDRGN_SRC/utils/analysis/analyze.sh [WORKDIR] [EPOCH] # Use 0-based index for the epoch number

The principle components of the learned heteorgeneity can be visualized with the `get_z_pcs.py` script to get the z-values along the PCS, followed by the `eval_decoder.py` script to generate the volumes.

    $ python $CDRGN_SRC/utils/analysis/get_z_pcs.py -h

## Fully unsupervised reconstruction

Please reach out to Ellen Zhong (zhonge[at]mit[dot]edu) if you'd like to collaborate on reconstruction of highly heterogeneous datasets where a consensus reconstruction is unavailable.  

## Contact

More documentation and tutorials to come! Bugs reports, feature requests, or general usage feedback to zhonge[at]mit[dot]edu.

