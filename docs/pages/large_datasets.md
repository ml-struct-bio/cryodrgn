# cryoDRGN for large datasets

CryoDRGN by default loads the entire dataset into memory for fast data access during training. However, large cryo-EM datasets can easily exceed the amount of memory available on standard workstations. For these datasets that do not fit into memory, `cryodrgn train_vae` can be run with the additional `--lazy` flag, which loads images on-the-fly instead of all at once at the beginning of training. This can, however, be slow due to the filesystem access pattern for on-the-fly image loading, especially if the data is not located on a SSD drive.

To improve processing speed of datasets when using the `--lazy` flag, you may want to bump up the `--max-threads` parameter to `8` or `16`, depending on the number of CPU cores available to you. This parameter will cause the input particle stack to be loaded using `max_threads` threads. In our tests, processing a large dataset on a spinning hard drive by specifying `--max-threads 16` improved the processing speed by a factor of 3.

```bash
# Parse pose information as usual, specifying the refinement box size with -D
cryodrgn parse_pose_csparc P10_J712_particles_exported.cs \
		-D 256 \
		-o data/pose.pkl

# Parse CTF information as usual
cryodrgn parse_ctf_csparc P10_J712_particles_exported.cs -o data/ctf.pkl

# Run cryoDRGN with extra flags --lazy and --max-threads
cryodrgn train_vae P10_J712_particles_exported.cs \
		--datadir P10/exports/groups/P10_J628_particles/J626/extract \
		--ctf data/ctf.pkl \
		--poses data/pose.pkl \
		--zdim 8 \
		-n 50 \
		--lazy \
		--max-threads 16 \
		-o 00_vae128 >> 00.log
```

## Numbers

Some numbers for training on a 1,375,854, 128x128 particle dataset (86 GB)

Baseline:

- 607 GB maximum memory requirement
- 18.5 min to load the dataset in `cryodrgn train_vae`

With `--lazy` and `--max-threads 16`:

- 200 GB maximum memory requirement
- 3.2 min to load the dataset in `cryodrgn train_vae`

On a single Nvidia V100 GPU, this dataset trained in approximately 2h,3min per epoch (large 1024x3 model) when fully loaded into memory. Training with on-the-fly data loading (`--lazy`) was 4x slower, though this can vary widely depending on your filesystem/network. Recent tests on cached filesystems do not have a large penalty for `--lazy` image loading.

## Still too large

If your dataset is still too large to load into memory, we recommend training on a subset of the images such that the dataset can fit into memory (e.g. split your dataset into two halves and run independent training jobs on each half). A random selection of a subset of your dataset can be generated with the utility `cryodrgn_utils select_random`:

```
# select 200k random particles out of a dataset containing 1,375,854 particles
(cryodrgn) $ cryodrgn_utils select_random 1375854 -n 200000 -o ind200k.pkl
```
