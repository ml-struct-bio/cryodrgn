# :snowflake::dragon: cryoDRGN: Deep Reconstructing Generative Networks for cryo-EM heterogeneous reconstruction

CryoDRGN is a neural network based algorithm for heterogeneous cryo-EM reconstruction. In particular, the method models a *continuous* distribution over 3D structures by using a neural network based representation for the volume.

## Manuscript:

CryoDRGN: Reconstruction of heterogeneous structures from cryo-electron micrographs using neural networks.
Ellen D. Zhong, Tristan Bepler, Bonnie Berger*, Joseph H. Davis*.
https://www.biorxiv.org/content/10.1101/2020.03.27.003871v1

Reconstructing continuous distributions of 3D protein structure from cryo-EM images.
Ellen D. Zhong, Tristan Bepler, Joseph H. Davis*, Bonnie Berger*.
ICLR 2020, https://arxiv.org/abs/1909.05215

## Installation/dependencies:

Until the cryoDRGN conda/pip package is available, for now, git clone the source code and install the following dependencies with anaconda, replacing the cudatoolkit version as necessary:

    # Create conda environment
    conda create --name cryodrgn python=3.7
    conda activate cryodrgn

    # Install dependencies
    conda install pytorch cudatoolkit=10.1 -c pytorch # Replace cudatoolkit version if needed
    conda install pandas

    # Clone source code and install
    git clone https://github.com/zhonge/cryodrgn.git
    cd cryodrgn
    git checkout 0.2.0
    python setup.py install --user 

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
    
    cryodrgn downsample [input particle stack] -D 128 -o particles.128.mrcs

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

Test that pose and CTF parameters were parsed correctly using the voxel-based backprojection script.
The goal is to quickly verify that there are no major problems with the extracted values and that the output structure resembles the structure from the consensus reconstruction before beginning training.

Example usage: 

    $ cryodrgn backproject_voxel projections.128.mrcs \
            --poses pose.pkl \
            --ctf ctf.pkl \
            --invert-data \ # Invert sign of dataset; Use if particles are white on black
            -o backproject.128.mrc

Check that the output structure `backproject.128.mrc` resembles the structure from the consensus reconstruction. 
It will not match exactly as the `backproject_voxel` command backprojects phase-flipped particles onto the voxel grid, and by default only uses the first 10k images. If the structure is too noisy, you can increase the number of images that are used with the `--first` argument.

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

Example command to train a 1-D latent variable cryoDRGN model for 50 epochs on an image dataset `projections.128.mrcs` with poses `pose.pkl` and ctf parameters `ctf.pkl`:

    $ cryodrgn train_vae projections.128.mrcs 
            --poses pose.pkl \
            --ctf ctf.pkl \
            --zdim 1 -n 50 \
            -o 00_vae128_z1
            
* Results will be saved in the specified directory `00_vae128_z10`.
            
Example command to train a 10-D latent variable cryoDRGN model for 50 epochs on an image dataset `projections.128.mrcs` with poses `pose.pkl` and ctf parameters `ctf.pkl`:

    $ cryodrgn train_vae projections.128.mrcs 
            --poses pose.pkl \
            --ctf ctf.pkl \
            --zdim 10 -n 50 \
            -o 01_vae128_z10

Example command to train a 10-D latent variable cryoDRGN model for 20 epochs on an image dataset `projections.256.mrcs` with poses `pose.pkl` and ctf parameters `ctf.pkl`:

    $ cryodrgn train_vae projections.256.mrcs 
            --poses pose.pkl \
            --ctf ctf.pkl \
            --zdim 10 -n 20 \
            --qdim 1024 --qlayers 3 --pdim 1024 --players 3 \
            -o 02_vae256_z10

The number of epochs `-n` refers to the number of full passes through the dataset for training, and should be modified depending on the number of particles in the dataset. For a 100k particle dataset, the above settings required ~10 min per epoch for D=128 images, and ~1 hour per epoch for D=256 images. 

If you would like to train longer, a training job can be extended with the `--load /path/to/weights.{last_epoch}.pkl` argument. For example to extend the training of the previous example to 50 epochs:

    $ cryodrgn train_vae projections.256.mrcs
            --poses pose.pkl \
            --ctf ctf.pkl \
            --zdim 10 -n 50 \
            --qdim 1024 --qlayers 3 --pdim 1024 --players 3 \
            -o 02_vae256_z10
            --load 02_vae_256_z10/weights.19.pkl # 0-based indexing

Note: while these settings worked well for the datasets we've tested, they are highly experimental for the general case as different datasets have diverse sources of heterogeneity. Please reach out to the authors with questions/consult -- we'd love to learn more.

### Local pose refinement -- BETA!

Depending on the quality of the consensus reconstruction, image poses may contain errors.
Image poses may be *locally* refined using the `--do-pose-sgd` flag. More details on this method to come!

### 6. Analysis of results

Once the model has finished training, the output directory will contain a configuration file `config.pkl`, neural network weights `weights.pkl`, image poses (if performing pose sgd) `pose.pkl`, and the predicted latent encoding for each image `z.pkl`. There is also a `run.log` file containing the average loss for each epoch which can be used to plot the learning curve for model training.




To analyze and visualize the learned latent space, use the scripts in the `utils/analysis` subdirectory. 

Additional structures may be generated using the `gen_volumes` script. For example:

    $ python $CDRGN_SRC/eval_decoder.py [WORKDIR]/weights.pkl --config [WORKDIR]/config.pkl -z ZVALUE -o reconstruct.sample1.mrc

Or to generate a trajectory using values of z given in a file `zvalues.txt`:

    $ python $CDRGN_SRC/eval_decoder.py [WORKDIR]/weights.pkl --config [WORKDIR]/config.pkl --zfile zvalues.txt -o [WORKDIR]/trajectory

For models with higher dimensional latent variables (zD>1), a shell script is provided which samples 20 structures from the latent space, runs UMAP dimensionality reduction, and creates a template jupyter-notebook in the working directory, which may be used for interactive visualization of the results:

    $ $CDRGN_SRC/utils/analysis/analyze.sh [WORKDIR] [EPOCH] # Use 0-based index for the epoch number

The principle components of the learned heteorgeneity can be visualized with the `get_z_pcs.py` script to get the z-values along the PCS, followed by the `eval_decoder.py` script to generate the volumes.

    $ python $CDRGN_SRC/utils/analysis/get_z_pcs.py -h

## Fully unsupervised reconstruction

Please reach out to Ellen Zhong (zhonge[at]mit[dot]edu) if you'd like to collaborate on reconstruction of highly heterogeneous datasets where a consensus reconstruction is unavailable.  

## Contact

More documentation and tutorials to come! Bugs reports, feature requests, or general usage feedback to zhonge[at]mit[dot]edu.

