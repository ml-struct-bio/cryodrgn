![pypi-downloads](https://img.shields.io/pypi/dm/cryodrgn?style=flat&label=PyPI%20Downloads&logo=pypi&logoColor=%233775A9&labelColor=%23FFF8EC)
![stable-release](https://img.shields.io/pypi/v/cryodrgn?style=flat&logo=pypi&logoColor=%233775A9&logoSize=auto&label=stable%20release&labelColor=%23FFF8EC)
![beta-release](https://img.shields.io/pypi/v/cryodrgn?pypiBaseUrl=https%3A%2F%2Ftest.pypi.org&style=flat&logo=pypi&logoColor=%233775A9&logoSize=auto&label=beta%20release&labelColor=%23FFF8EC)
![grading](https://img.shields.io/codefactor/grade/github/ml-struct-bio/cryodrgn/main?style=flat&logo=codefactor&logoColor=%23F44A6A&logoSize=auto&label=CodeFactor%20Grade&labelColor=%23FFF8EC)
![ci-test](https://github.com/ml-struct-bio/cryodrgn/actions/workflows/tests.yml/badge.svg)


# :snowflake::dragon: cryoDRGN: Deep Reconstructing Generative Networks for cryo-EM and cryo-ET heterogeneous reconstruction

CryoDRGN is a neural network based algorithm for heterogeneous cryo-EM reconstruction. In particular, the method models
a *continuous* distribution over 3D structures by using a neural network based representation for the volume.


## Documentation

The latest documentation for cryoDRGN is available in our [user guide](https://ez-lab.gitbook.io/cryodrgn/), including an overview and walkthrough of
cryoDRGN installation, training and analysis. A brief quick start is provided below.

For any feedback, questions, or bugs, please file a Github issue or start a Github discussion.


### New in Version 4.x

* [NEW] integration of CryoDRGN-AI ab-initio reconstruction method
* [NEW] the `cryodrgn train` command for access to all reconstruction methods, including fixed poses and
  ab-initio models, using configuration files for model parameters instead of command-line arguments


### New in Version 3.x

* [NEW] official release of [cryoDRGN-ET](https://www.biorxiv.org/content/10.1101/2023.08.18.553799v1) for heterogeneous subtomogram analysis.
* [NEW] Heterogeneous reconstruction of subtomograms. See documentation [on gitbook](https://ez-lab.gitbook.io/cryodrgn/)
* Updated `cryodrgn backproject_voxel` for voxel-based homogeneous reconstruction
* Major refactor of dataset loading for handling large datasets
* [NEW] `cryodrgn plot_classes` for analysis visualizations colored by a given set of class labels
* implementing [automatic mixed-precision training](https://pytorch.org/docs/stable/amp.html) for ab-initio
  reconstruction for 2-4x speedup
* support for RELION 3.1 .star files with separate optics tables, np.float16 number formats used in RELION .mrcs outputs
* `cryodrgn backproject_voxel` produces cryoSPARC-style FSC curve plots with phase-randomization correction of
  automatically generated tight masks
* official support for Python 3.11 and 3.12


## Previous versions

<details><summary>Version 3.3</summary><ul>
  <li>[NEW] <code>cryodrgn direct_traversal</code> to generate interpolations in the conformation latent space
between two points</li>
  <li>support for .txt files in <code>write_star</code></li>
  <li>adding <code>--datadir</code> to <code>cryodrgn abinit_homo</code> for use with .star files</li>
  <li>fixing various bugs in <code>backproject_voxel</code>, Jupyter demonstration notebooks</li>
  <li>support for TestPyPI beta release deployments via <code>pip</code></li>
</ul></details>

<details><summary>Version 3.2</summary><ul>
  <li>[NEW] <code>cryodrgn_utils clean</code> for removing extraneous output files from completed experiments</li>
  <li>[NEW] <code>cryodrgn_utils fsc</code>, <code>cryodrgn_utils fsc_plot</code>,
<code>cryodrgn_utils gen_mask</code>adapted from existing scripts — for calculating FSCs, plotting them, and
generating masks for volumes respectively</li>
  <li><code>cryodrgn backproject_voxel</code> now produces half-maps and a half-map FSC</li>
  <li>fixing <code>filter_star</code> to accept tilt series as well</li>
  <li>fixing assorted bugs in e.g. <code>write_star</code>, <code>invert_constrast</code>, and
<code>train_vae</code> (see release notes)
</ul></details>

<details><summary>Version 3.1</summary><ul>
  <li><code>cryodrgn filter</code> interface for interactive filtering of particles as an alternative to the
cryoDRGN_filter Jupyter notebook</li>
</ul></details>

<details><summary>Version 2.3</summary><ul>
  <li>Model configuration files are now saved as human-readable config.yaml files
(https://github.com/zhonge/cryodrgn/issues/235)</li>
  <li>Fix machine stamp in output .mrc files for better compatibility with downstream tools
(https://github.com/zhonge/cryodrgn/pull/260)</li>
  <li>Better documentation of help flags in ab initio reconstruction tools
(https://github.com/zhonge/cryodrgn/issues/258)</li>
  <li>[FIX] By default, window images in <code>cryodrgn abinit_homo</code> (now consistent with other reconstruction
tools) (https://github.com/zhonge/cryodrgn/issues/258)</li>
  <li>[FIX] Reduce memory usage when using <code>--preprocessed</code>
and <code>--ind</code> (https://github.com/zhonge/cryodrgn/pull/272)</li>
</ul></details>

<details><summary>Version 2.2</summary><ul>
  <li>[NEW] Tools for ab initio homogeneous and heterogeneous reconstruction:

    (cryodrgn) $ cryodrgn abinit_homo -h
    (cryodrgn) $ cryodrgn abinit_het -h
  </li>

  <li>[NEW] Utils function for writing cryoSPARC .cs/.csg files
[to reimport data into cryoSPARC](https://github.com/zhonge/cryodrgn/issues/150#issuecomment-1465094751):


    (cryodrgn) $ cryodrgn_utils write_cs
  </li>
  <li>[Improved plotting](https://github.com/zhonge/cryodrgn/issues/219) in <code>cryodrgn analyze</code></li>
  <li>Many codebase improvements with open-source software development practices (e.g. continuous integration tests,
black, flake8, pyright, logging, and PyPi packaging).</li>
</ul></details>

<details><summary>Version 1.1.x</summary>
Updated default parameters for <code>cryodrgn train_vae</code> with modified positional encoding, larger model
architecture, and accelerated mixed-precision training turned on by default:
<ul>
  <li>Mixed precision training is now turned on by default (Use <code>--no-amp</code> to revert to single precision
training)</li>
  <li>Encoder/decoder architecture is now 1024x3 by default (Use <code>--enc-dim 256</code> and <code>--dec-dim
256</code> to revert)
</li>
  <li>Gaussian Fourier featurization for faster training and higher resolution density maps (Use <code>--pe-type
geom_lowf</code>
to revert)</li>
</ul></details>

<details><summary>Version 1.0.x</summary>
The official version 1.0 release. This version introduces several new tools for analysis of the reconstructed ensembles,
and adds functionality for calling utility scripts with <code>cryodrgn_utils {command}</code>.
<ul>
  <li>NEW: <code>cryodrgn analyze_landscape</code> and <code>cryodrgn analyze_landscape_full</code> for automatic
assignment of classes and conformational landscape visualization. Documentation for this new feature is here:
https://www.notion.so/cryodrgn-conformational-landscape-analysis-a5af129288d54d1aa95388bdac48235a.</li>
  <li>NEW: Faster training and higher resolution model
with Gaussian Fourier featurization (Use <code>--pe-type gaussian</code>)</li>
  <li>NEW: <code>cryodrgn_utils {command} -h</code> for standalone utility scripts</li>
  <li>NEW: <code>cryodrgn_utils write_star</code> for converting cryoDRGN particle selections to .star files</li>
  <li>Add pytorch native mixed precision training and fix support for pytorch 1.9+</li>
</ul></details>

<details><summary>Version 0.3.4</summary><ul>
    <li>FIX: Bug in write_starfile.py when provided particle stack is chunked (.txt file)</li>
    <li>Support micrograph coordinates and additional column headers to write_starfile.py</li>
    <li>New helper scripts: analyze_convergence.py (<i>in beta testing</i>) contributed by <a href="bmp@mit.
edu">Barrett
Powell</a> (thanks!) and make_random_selection.py for splitting up particle stacks for training</li>
</ul></details>

<details><summary>Version 0.3.3</summary><ul>
    <li>Faster image preprocessing and smaller memory footprint</li>
    <li>New: <code>cryodrgn preprocess</code> for large datasets (<i>in beta testing</i> - see
<a href="https://www.notion.so/cryodrgn-preprocess-d84a9d9df8634a6a8bfd32d6b5e737ef">this Notion doc</a> for details)
</li>
    <li>* Known <a href="https://github.com/zhonge/cryodrgn/issues/66">issue</a> with PyTorch version 1.9+</li>
</ul></details>

<details><summary>Version 0.3.2</summary><ul>
    <li>New: cryoDRGN_filtering.ipynb for interactive filtering and selection of images from the dataset</li>
    <li>New: <code>cryodrgn view_config</code></li>
    <li>Minor performance improvements and compatibility fixes</li>
</ul></details>

<details><summary>Version 0.3.1</summary><ul>
    <li>New: Script write_starfile.py to convert (filtered) particle selection to a .star file</li>
    <li>More visualizations in <code>cryodrgn analyze</code></li>
</ul></details>

<details><summary>Version 0.3.0</summary><ul>
    <li>New: GPU parallelization with flag <code>--multigpu</code></li>
    <li>New: Mode for accelerated mixed precision training with flag <code>--amp</code>, available for NVIDIA
tensor core GPUs</li>
    <li>Interface update: renamed encoder arguments <code>--qdim</code> and <code>--qlayers</code> to
<code>--enc-dim</code> and <code>--enc-layers</code>; renamed decoder arguments <code>--pdim</code> and
<code>--players</code> to <code>--dec-dim</code> and <code>--dec-layers</code></li>
    <li>Argument default changes: flipped the default for <code>--invert-data</code> to True by default, and
flipped the default for <code>--window</code> to True by default</li>
    <li>Updated training recommendations in below quick start guide</li>
    <li>Updates to cryodrgn analyze: more visualizations, ordering kmeans volumes according to distances in latent
space (previously random), and more features for particle selection and filtering in the Jupiter notebook</li>
</ul></details>

<details><summary>Version 0.2.1</summary><ul>
    <li>New: Parsing of RELION 3.1 files</li>
    <li>Fix: Compatibility with pytorch 1.5</li>
</ul></details>

<details><summary>Version 0.2.0</summary><ul>
    <li>New interface and proper python packaging with setup.py. This version has identical functionality and
argument usage as previous versions, however tools are now available from a common entry point. See:

```$ cryodrgn -h```

</li>
<li>New analysis pipeline <code>cryodrgn analyze</code></li>
<li>New latent space traversal scripts with
<code>cryodrgn graph_traversal</code> and <code>cryodrgn pc_traversal</code>.</li>
</ul></details>


## Installation

The alpha releases of `cryodrgn v4.x` may be installed via `pip`, and we recommend installing `cryodrgn` in a clean
conda environment.

    # Create and activate conda environment
    (base) $ conda create --name cryodrgn python=3.10
    (cryodrgn) $ conda activate cryodrgn

    # install cryodrgn
    (cryodrgn) $ pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ \
                             --pre 'cryodrgn>4'

You can alternatively install a more stable version of `cryodrgn` using our primary release channel:

    (cryodrgn) $ pip install cryodrgn

More installation instructions are found in the [documentation](https://ez-lab.gitbook.io/cryodrgn/installation).

## Quickstart: heterogeneous reconstruction with consensus poses

### 1. Preprocess image stack

First resize your particle images using the `cryodrgn downsample` command.
We recommend first downsampling images to 128x128 since larger images can take much longer to train:

    $ cryodrgn downsample [input particle stack] -D 128 -o particles.128.mrcs

If there are memory issues with downsampling large particle stacks, add the `--chunk 10000` argument to
save images as separate `.mrcs` files of 10k images.

### 2. Parse image poses from a consensus homogeneous reconstruction

CryoDRGN expects image poses to be stored in a binary pickle format (`.pkl`). Use the `parse_pose_star` or
`parse_pose_csparc` command to extract the poses from a `.star` file or a `.cs` file, respectively.

Example usage to parse image poses from a RELION 3.1 starfile:

    $ cryodrgn parse_pose_star particles.star -o pose.pkl

Example usage to parse image poses from a cryoSPARC homogeneous refinement particles.cs file:

    $ cryodrgn parse_pose_csparc cryosparc_P27_J3_005_particles.cs -o pose.pkl -D 300

**Note:** The `-D` argument should be the box size of the consensus refinement (and not the downsampled
images from step 1) so that the units for translation shifts are parsed correctly.


### 3. Parse CTF parameters from a .star/.cs file

CryoDRGN expects CTF parameters to be stored in a binary pickle format (`.pkl`).
Use the `parse_ctf_star` or `parse_ctf_csparc` command to extract the relevant CTF parameters from a `.star` file
or a `.cs` file, respectively.

Example usage for a .star file:

    $ cryodrgn parse_ctf_star particles.star -o ctf.pkl

If the box size and Angstrom/pixel values are not included in the .star file under fields `_rlnImageSize` and
`_rlnImagePixelSize` respectively, the `-D` and `--Apix` arguments to `parse_ctf_star` should be used instead to
provide the original parameters of the input file (before any downsampling):

    $ cryodrgn parse_ctf_star particles.star -D 300 --Apix 1.03 -o ctf.pkl


### 4. (Optional) Test pose/CTF parameters parsing

Next, test that pose and CTF parameters were parsed correctly using the voxel-based backprojection script.
The goal is to quickly verify that there are no major problems with the extracted values and that the output structure
resembles the structure from the consensus reconstruction before training.

Example usage:

    $ cryodrgn backproject_voxel projections.128.mrcs --poses pose.pkl --ctf ctf.pkl -o backproject.128 --first 10000

The output structure `backproject.128/backproject.mrc` will not match the consensus reconstruction exactly
as the `backproject_voxel` command backprojects phase-flipped particles onto the voxel grid, and because here we
performed backprojection using only the first 10k images in the stack for quicker results.
If the structure is too noisy, we can try using more images with `--first` or the
entire stack instead (without `--first`).

**Note:** If the volume does not resemble your structure, you may need to use the flag `--uninvert-data`.
This flips the data sign (e.g. light-on-dark or dark-on-light), which may be needed depending on the
convention used in upstream processing tools.

### 5. Running cryoDRGN heterogeneous reconstruction

When the input images (.mrcs), poses (.pkl), and CTF parameters (.pkl) have been prepared, a cryoDRGN model
can be trained with the `cryodrgn train_vae` command:

    $ cryodrgn train_vae particles_128.mrcs --ctf ctf.pkl --poses pose.pkl -o 001_train-vae.128 \
                                            --zdim 8 -n 50 --dec-dim=128 --enc-dim=128

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

After validation, pose optimization, and any necessary particle filtering,
you can then train on the full resolution images!


## 6. Analysis of results

Once the model has finished training, the output directory will contain a configuration file `config.yaml`,
neural network weights `weights.pkl`, image poses (if performing pose sgd) `pose.pkl`,
and the latent embeddings for each image `z.pkl`.
The latent embeddings are provided in the same order as the input particles.
To analyze these results, use the `cryodrgn analyze` command to visualize the latent space and generate structures.
`cryodrgn analyze` will also provide a template jupyter notebook for further interactive visualization and analysis.

Example usage to analyze results from the direction `01_cryodrgn256` containing results after 25 epochs of training:

    $ cryodrgn analyze 01_cryodrgn256 24 --Apix 1.31 # 24 for 0-based indexing of epoch numbers


## References:

For a complete description of the method, see:

* CryoDRGN: reconstruction of heterogeneous cryo-EM structures using neural networks
Ellen D. Zhong, Tristan Bepler, Bonnie Berger*, Joseph H Davis*
Nature Methods 2021, https://doi.org/10.1038/s41592-020-01049-4 [pdf](https://ezlab.princeton.edu/assets/pdf/2021_cryodrgn_nature_methods.pdf)

CryoDRGN2's ab initio reconstruction algorithms were published at ICCV:

* CryoDRGN2: Ab Initio Neural Reconstruction of 3D Protein Structures From Real Cryo-EM Images
Ellen D. Zhong, Adam Lerer, Joseph H Davis, and Bonnie Berger
International Conference on Computer Vision (ICCV) 2021, [paper](https://openaccess.thecvf.com/content/ICCV2021/papers/Zhong_CryoDRGN2_Ab_Initio_Neural_Reconstruction_of_3D_Protein_Structures_From_ICCV_2021_paper.pdf)


## Contact

Please submit any bug reports, feature requests, or general usage feedback as a github issue or discussion.
