# ❄️🐉 cryoDRGN Introduction

This document contains a guide for installing and running **cryoDRGN** 🐉 ❄️. In particular, we follow the processing steps for **particle filtering** and **heterogeneous reconstruction** of the **assembling ribosome dataset (EMPIAR-10076)** used in [Zhong et al](https://www.nature.com/articles/s41592-020-01049-4). This is meant as a general guide — submission commands may need to be updated depending on your workstation or cluster setup.

For any feedback, issues, or typos, please file a Github [issue](https://github.com/zhonge/cryodrgn/issues) or send an email to the cryodrgn user [google group](https://groups.google.com/g/cryodrgn).

---

## Background

CryoDRGN is a neural network-based method for heterogeneous reconstruction. Instead of *discrete* methods like 3D classification that produce an ensemble of K density maps, cryoDRGN performs heterogeneous reconstruction by learning a *continuous distribution* of density maps parameterized by a coordinate-based neural network.

<iframe src="https://widgets.figshare.com/articles/21170578/embed?show_title=1" width="568" height="351" allowfullscreen frameborder="0"></iframe>

*Principal component trajectories and graph traversal trajectories of the pre-catalyic spliceosome. SI Video 4 from [Zhong et al 2021](https://www.nature.com/articles/s41592-020-01049-4)*

The inputs to a cryoDRGN training run are **1) extracted particle images**, **2) the CTF parameters** associated with each particle, and **3) poses** for each particle from a 3D refinement. Note that cryoDRGN treats the reconstruction as C1 (asymmetric). For a few thoughts on (pseudo-)symmetric complexes, see this [note](https://github.com/zhonge/cryodrgn/issues/21).

The final result of the software will be **1) latent embeddings** for each particle image in the form of a real-valued vector (usually denoted with z, and output as a `z.pkl` file by the software), and **2) neural network weights** modeling the distribution of density maps (parameterizing the function from z→V). Once trained, the software can reconstruct a 3D density map given a value of z.

How do you interpret the resulting distribution of structures? Since different datasets have diverse sources of heterogeneity (e.g. discrete vs. continuous), cryoDRGN contains a variety of automated and interactive tools to analyze the reconstructed distribution of structures. The starting point for analysis is the `cryodrgn analyze` pipeline, which generates a sample of 3D density maps and visualizations of the latent space. Specifically, the `cryodrgn analyze` pipeline will produce **1) N density maps** sampled from different regions of the latent space (N=20, by default), **2) continuous trajectories** along the principal components axes of the latent space embeddings, and **3) visualizations of the latent space** with PCA and UMAP**.**

CryoDRGN also provides interactive tools to further explore the learned ensemble, implemented as **Jupyter notebooks** with interactive widgets for visualizing the dataset, extracting particles, and generating more volumes. Additional tools are also available that can generate trajectories given user-defined end points and convert particle selections to `.star` files for further refinement in other tools. An overview of these functionalities will be demonstrated in the tutorial.

Furthermore, because the model is trained to reconstruct *image heterogeneity,* any non-structural image heterogeneity that is not captured by the image formation model **(e.g. junk particles and artifacts)** can be reflected in the latent embeddings. In practice, junk particles are often easily identified in the latent embeddings and can then be filtered out. A jupyter notebook is provided to filter particle stacks.

What settings should I use for training cryoDRGN networks? Common hyperparameters when training a cryoDRGN model are: **1) the size of the neural network**, which controls the capacity of the model, **2) the input image size**, which bounds the resolution information and greatly impacts the training speed and **3) the latent variable dimension**, which is the bottleneck layer that bounds the expressiveness of the model. The three parameters together all affect the expressiveness/complexity of the learned model. After exploring many real datasets, we provide reasonable defaults and recommended settings of these parameters for training.

## Input data requirements

- Extracted single particle images (in .mrcs/.cs/.star format), clean from edge, ice, or hot pixel artifacts
- A C1 consensus reconstruction with:
    - High quality CTF parameters
    - High quality image poses (particle alignments)

## Tutorial overview

See [cryoDRGN EMPIAR-10076 tutorial](empiar_tutorial.md) for a step-by-step guide for running cryoDRGN.

This walkthrough of cryoDRGN analysis of the **assembling ribosome dataset (EMPIAR-10076)** covers:

1. preprocessing of inputs,
2. initial cryoDRGN training and explanation of outputs,
3. particle filtering to remove junk particles,
4. high-resolution cryoDRGN training,
5. extracting particle subsets for traditional refinement, and
6. generation of trajectories.

For an abbreviated overview of the steps for running cryoDRGN, see the github [README](https://github.com/zhonge/cryodrgn)

A protocols paper that describes the analysis of the assembling ribosome dataset is now published. See [Kinman*, Powell*, Zhong* et al.](https://www.nature.com/articles/s41596-022-00763-x)

<iframe src="https://widgets.figshare.com/articles/21170908/embed?show_title=1" width="568" height="351" allowfullscreen frameborder="0"></iframe>

*SI Video 3 from [Zhong et al 2021](https://www.nature.com/articles/s41592-020-01049-4)*

## References

For a complete description of the method, see our paper here:

**CryoDRGN: reconstruction of heterogeneous cryo-EM structures using neural networks**

Ellen Zhong, Tristan Bepler, Bonnie Berger*, Joey Davis*

Nature Methods 2021, [https://doi.org/10.1038/s41592-020-01049-4](https://doi.org/10.1038/s41592-020-01049-4)

---

An earlier version of this work appeared at the International Conference of Learning Representations (ICLR):

**Reconstructing continuous distributions of protein structure from cryo-EM images**

Ellen Zhong, Tristan Bepler, Joey Davis*, Bonnie Berger*

ICLR 2020, Spotlight, [https://arxiv.org/abs/1909.05215](https://arxiv.org/abs/1909.05215)
