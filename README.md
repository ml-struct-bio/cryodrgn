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

Until the cryoDRGN conda package is available, for now, git clone the source code and install the following dependencies with anaconda, replacing the cudatoolkit version as necessary:

    conda create --name cryodrgn python=3.7
    conda activate cryodrgn
    conda install pytorch torchvision cudatoolkit=10.1 -c pytorch
    conda install pandas

Additional requirements for latent space analysis and interactive visualization:

    export CDRGN_SRC="path/to/git/repo"
    conda install seaborn scikit-learn 
    conda install -c conda-forge umap-learn
    conda install -c conda-forge jupyterlab
    pip install ipywidgets
    jupyter nbextension enable --py widgetsnbextension
    pip install cufflinks

## Quickstart: heterogeneous reconstruction with consensus poses

### 1. Preprocess image stack

Training cryoDRGN networks has not been tested on image sizes above D=256. If your images are larger than D=256, use the following utility to downsample the images:

    $ python $CDRGN_SRC/utils/fouriershrink.py [input particle stack] -D 256 -o [output .mrcs]

It is also recommended to create image stacks at lower resolution (e.g. D=128) for initial testing and pilot experiments with cryoDRGN.

    $ python $CDRGN_SRC/utils/fouriershrink.py [input particle stack] -D 128 -o [output .mrcs]
    
The input file format can be a single `.mrcs` file, a `.txt` file containing paths to multiple `.mrcs` files, a `.star` file, or a cryoSPARC `.cs` file. For the latter two options, if the relative paths to the `.mrcs` are broken, the argument `--datadir` can be used to supply the path to where the `.mrcs` files are located. 

If there are memory issues with large particle stacks, add the `--chunk 10000` argument to save out images as separate `.mrcs` files of 10k images. 

### 2. Parse image poses from a consensus homogeneous reconstruction

To parse image poses from a RELION starfile:
    
    $ python $CDRGN_SRC/utils/parse_pose_star.py particles.star -o pose.pkl -D 300

To parse image poses from a cryoSPARC homogeneous refinement particles.cs file:

    $ python $CDRGN_SRC/utils/parse_pose_csparc.py cryosparc_P27_J3_005_particles.cs -o pose.pkl -D 300

The `-D` argument should be set to the box size of the original reconstruction (before any downsampling). 

### 3. Parse CTF parameters from a .star/.cs file

CryoDRGN currently loads CTF parameters from a binary pickle format (`.pkl`). Use the utility script `parse_ctf_star.py` or `parse_ctf_csparc.py` to extract the relevant CTF parameters from a `.star` file or a `.cs` file, respectively.

Example usage:
    
    # .star file
    $ python $CDRGN_SRC/utils/parse_ctf_star.py particles.star -D 300 --Apix 1.03 -o ctf.pkl
    # .cs file
    $ python $CDRGN_SRC/utils/parse_ctf_csparc.py cryosparc_P27_J3_005_particles.cs -o ctf.pkl

The `-D` and `--Apix` arguments should be set to the box size and Angstrom/pixel of the original `.mrcs` file (before any downsampling). 


### 4. Test pose/CTF parameters were parsed correctly

Test that pose and CTF parameters were parsed correctly using the voxel-based backprojection script:

    $ python $CDRGN_SRC/backproject_voxel.py projections.128.mrcs \
            --poses pose.pkl \
            --ctf ctf.pkl \
            --invert-data \ # Invert sign of dataset; Use if particles are white on black
            -o backproject.128.mrc

Check that the output structure `backproject.128.mrc` resembles the structure from the consensus reconstruction. 
It will not match exactly as the `backproject_voxel.py` script backprojects phase-flipped particles onto the voxel grid, and by default only uses the first 10k images. Increase the number of images that are used with the `--first` argument.

### 5. Running cryoDRGN heterogeneous reconstruction

When the input image stack (.mrcs), image poses (.pkl), and CTF parameters (.pkl) have been prepared, a cryoDRGN model can be trained with following script:

    $ python $CDRGN_SRC/vae_het.py -h
    
Many of the parameters of this script have sensisible defaults. The required arguments are:

* an input image stack (`.mrcs`)
* `--poses`, image poses (`.pkl`)
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

    $ python $CDRGN_SRC/vae_het.py projections.128.mrcs 
            --poses pose.pkl \
            --ctf ctf.pkl \
            --zdim 1 -n 50 \
            -o 00_vae128_z1
            
* Results will be saved in the specified directory `00_vae128_z10`.
            
Example command to train a 10-D latent variable cryoDRGN model for 50 epochs on an image dataset `projections.128.mrcs` with poses `pose.pkl` and ctf parameters `ctf.pkl`:

    $ python $CDRGN_SRC/vae_het.py projections.128.mrcs 
            --poses pose.pkl \
            --ctf ctf.pkl \
            --zdim 10 -n 50 \
            -o 01_vae128_z10

Example command to train a 10-D latent variable cryoDRGN model for 20 epochs on an image dataset `projections.256.mrcs` with poses `pose.pkl` and ctf parameters `ctf.pkl`:

    $ python $CDRGN_SRC/vae_het.py projections.256.mrcs 
            --poses pose.pkl \
            --ctf ctf.pkl \
            --zdim 10 -n 20 \
            --qdim 1024 --qlayers 3 --pdim 1024 --players 3 \
            -o 02_vae256_z10

The number of epochs `-n` should be modified depending on the number of particles in the dataset. The above parameters led to reasonable training times on datasets with ~100-300k particles (~hours on D=128 images, ~days on D=256 images). 

Note: while these settings worked well for the datasets we've tested, they are highly experimental for the general case as different datasets have diverse sources of heterogeneity. Please reach out to the authors with questions/consult -- we'd love to learn more.

### Local pose refinement -- BETA!

Depending on the quality of the consensus reconstruction, image poses may contain errors.
Image poses may be *locally* refined using the `--do-pose-sgd` flag. 

### 6. Analysis of results

Once the model has finished training, the output directory will contain a configuration file `config.pkl`, neural network weights `weights.pkl`, image poses (if performing pose sgd) `pose.pkl`, and the predicted latent encoding for each image `z.pkl`. 

To analyze and visualize the learned latent space, use the scripts in the `utils/analysis` subdirectory. 

Additional structures may be generated using the `eval_decoder.py` script. For example:

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

