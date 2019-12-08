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
    conda install seaborn scikit-learn 
    conda install -c conda-forge umap-learn
    pip install --user healpy

## Quickstart: heterogeneous reconstruction with consensus alignments

### 1. Preprocess image stack

Training cryoDRGN networks has not been tested on image sizes above D=256. If your images are larger than D=256, use the following utility to downsample the images:

    $ SRC=[path to source code]
    $ python $SRC/utils/fouriershrink.py [input particle stack] -D 256 -o [output particle stack] --out-png projections.256.png

It is also recommended to create image stacks at lower resolution (e.g. D=80, 160) for initial testing and pilot experiments with cryoDRGN.

    $ python $SRC/utils/fouriershrink.py [input particle stack] -D 160 -o [output particle stack] --out-png projections.160.png
    $ python $SRC/utils/fouriershrink.py [input particle stack] -D 80 -o [output particle stack] --out-png projections.80.png

### 2. Parse alignments from a consensus homogeneous reconstruction

* Parse alignments from RELION starfile

    $ python $SRC/utils/parse_star_alignments.py particles.star -o consensus

* Parse alignments from cryoSPARC particles.cs file

    $ python $SRC/utils/parse_cs_alignments.py cryosparc_P27_J3_005_particles.cs --homorefine -o consensus

* Test that alignments were parsed correctly using the voxel-based backprojection script:

    $ python $SRC/backproject_voxel.py -h

### 3. Parse CTF parameters from a .star file

CryoDRGN currently takes CTF parameters in a pickle format (.pkl). Use the utility script `parse_ctf.py` to extract the relevant CTF parameters from a .star file. 

Example usage:

    $ python $SRC/utils/parse_ctf.py particles.star -N 101845 -D 256 --Apix 1.7 -o ctf.256.pkl

Note: the pixel size is saved in the CTF pickle. Run this script multiple times for each pixel size if cryoDRGN will be run on variable image size particles.

    $ python $SRC/utils/parse_ctf.py particles.star -N 101845 -D 160 --Apix 2.72 -o ctf.160.pkl
    $ python $SRC/utils/parse_ctf.py particles.star -N 101845 -D 80 --Apix 5.44 -o ctf.80.pkl

### 4. Running cryoDRGN heterogeneous reconstruction

When the input image stack (.mrcs), image poses (.pkl), and CTF parameters (.pkl) have been prepared, the cryoDRGN networks can be trained with following script:

    $ python $SRC/vae_het.py -h
    usage: vae_het.py [-h] -o OUTDIR --zdim ZDIM --poses [POSES [POSES ...]]
                      [--tscale TSCALE] [--ctf pkl] [--load LOAD]
                      [--checkpoint CHECKPOINT] [--log-interval LOG_INTERVAL] [-v]
                      [--seed SEED] [--invert-data] [--window] [--ind IND]
                      [--tilt TILT] [--tilt-deg TILT_DEG] [-n NUM_EPOCHS]
                      [-b BATCH_SIZE] [--wd WD] [--lr LR] [--beta BETA]
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
      particles             Particles (.mrcs)
    
    optional arguments:
      -h, --help            show this help message and exit
      -o OUTDIR, --outdir OUTDIR
                            Output directory to save model
      --zdim ZDIM           Dimension of latent variable
      --poses [POSES [POSES ...]]
                            Image rotations and optionally translations (.pkl)
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
      --invert-data         Invert data sign
      --window              Real space windowing of dataset
      --ind IND             Filter particle stack by these indices
    
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
      --pose-lr POSE_LR     Learning rate for pose optimizer (default: 0.0001)
    
    Encoder Network:
      --qlayers QLAYERS     Number of hidden layers (default: 10)
      --qdim QDIM           Number of nodes in hidden layers (default: 128)
      --encode-mode {conv,resid,mlp,tilt}
                            Type of encoder network (default: resid)
      --enc-mask ENC_MASK   Circular mask of image for encoder (default: D/2; -1
                            for no mask)
      --use-real            Use real space image for encoder (for convolutional
                            encoder)
    
    Decoder Network:
      --players PLAYERS     Number of hidden layers (default: 10)
      --pdim PDIM           Number of nodes in hidden layers (default: 128)
      --pe-type {geom_ft,geom_full,geom_lowf,geom_nohighf,linear_lowf,none}
                            Type of positional encoding (default: geom_lowf)
      --domain {hartley,fourier}
                            Decoder representation domain (default: fourier)

### Example usage:

Example command to train a 10-D latent variable cryoDRGN model for 20 epochs on an image dataset `projections.256.mrcs` with poses `consensus.rot.pkl, consensus.trans.pkl` and ctf parameters `ctf.256.pkl`:

    $ python $SRC/vae_het.py projections.256.mrcs \
            --poses consensus.rot.pkl consensus.trans.pkl \
            --tscale .8 \ 
            --ctf ctf.256.pkl \
            --zdim 10 \
            --do-pose-sgd \ 
            -n 20 \
            --pdim 1000 --players 3 \
            --qdim 1000 --qlayers 3 \
            -o 00_vae256_z10

* The encoder and decoder networks in this model contain 3 layers of dimension 1000. 
* Using the `--do-pose-sgd` flag will lead to *local* refinement of poses.
* Results will be saved in the specified directory `00_vae256_z10`.
* NOTE: Since translations are provided in units of pixels, `--tscale` may be used to renormalize translations e.g. if the consensus poses were obtained with a different box size.

### Recommended settings:

Pilot experiments:
* Run with small model, small image sizes, and zD=1
* Run with small model, small image sizes, and zD=10

After validation, pose optimization, and any necessary particle filtering, train the full resolution image stack (up to D=256) with a large model:
* Run large model, large image sizes, zD=10

Note these settings are highly experimental.

### Analysis

Once the model has finished training, the working directory will contain a configuration file `config.pkl`, neural network weights `weights.pkl`, image poses (if performing pose sgd) `pose.pkl`, and the predicted latent encoding for each image `z.pkl`. The directory will also contain a volume `reconstruct.mrc` evaluated at the mean z value of the dataset. 

To analyze and visualize the learned latent space, use the scripts in the `utils/analysis` subdirectory (TODO: Document). Additional structures may be generated using the `eval_decoder.py` script. For example:

    $ python $SRC/eval_decoder.py [WORKDIR]/weights.pkl --config [WORKDIR]/config.pkl -z ZVALUE -o reconstruct.sample1.mrc

Or to generate a trajectory using values of z given in a file `zvalues.txt`:

    $ python $SRC/eval_decoder.py [WORKDIR]/weights.pkl --config [WORKDIR]/config.pkl --zfile zvalues.txt -o [WORKDIR]/trajectory

## Fully unsupervised reconstruction

Please reach out to Ellen Zhong (zhonge[at]mit[dot]edu) if you'd like to collaborate on reconstruction of highly heterogeneous datasets where a consensus reconstruction is unavailable.  

## Contact

More documentation and tutorials to come! Bugs reports, feature requests, or general usage feedback to zhonge[at]mit[dot]edu.

