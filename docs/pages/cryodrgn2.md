# CryoDRGN2 quickstart

There are two commands for ab initio reconstruction, `cryodrgn abinit_homo` and `cryodrgn abinit_het` for homogeneous and heterogeneous ab initio reconstruction, respectively:

```bash
# homogeneous ab initio reconstruction
(cryodrgn) $ cryodrgn abinit_homo -h

# heterogeneous ab initio reconstruction
(cryodrgn) $ cryodrgn abinit_het -h
```

### Setup

- Downsample your particles to a box size of 128 either with `cryodrgn downsample` or with other tools.
- If you have a large dataset (>500k images), we recommend training on a subset of particles for initial testing. Use `cryodrgn_utils select_random` to select a random subset of particles.

    ```bash
    # get a random selection of 200k particles from a dataset of 1,423,124 particles
    (cryodrgn) $ cryodrgn_utils select_random 1423124 -n 200000 -o ind200k.pkl
    ```

    - You can then train on only the random subset with the argument `--ind ind200k.pkl`
- For reference, ab initio heterogeneous reconstruction on a dataset containing 218k 128x128 particles took 20 hours to train on a single V100 GPU.

### Example usage

```bash
# homogeneous reconstruction
(cryodrgn) $ cryodrgn abinit_homo [particles] --ctf [ctf.pkl] -o [output_directory]  >> output.log

# heterogeneous reconstruction
(cryodrgn) $ cryodrgn abinit_het [particles] --ctf [ctf.pkl] --zdim 8 -o [output_directory]  >> output.log
```

### Note on training settings

- The default translational search extent is +/- 10 pixels (`--t-extent 10`). If your particles are not well-centered, you can use a wider search extent, e.g. +/- 40 pixels ( `--t-extent 40`).
- Poses are updated every 5 epochs (`--ps-freq 5`) to alternate between pose search (slow) and standard cryodrgn1 training (fast) using the last iteration's poses.
- The default pose search settings are not tuned for high accuracy alignments (a tradeoff of accuracy vs. compute speed). You can increase the resolution of the pose search with `
- The default training time is 30 epochs. A typical use case is to run for 30 epochs, check the results (`cryodrgn analyze`), then extend training to 60 epochs. You can extend by rerunning with `-n 60 --load latest`. If your dataset is very large, you may want to reduce the pose search freqency `--ps-freq` and the number of epochs `-n`.
- During training, pose search epochs will get successively slower. This is because the parameter `--l-ramp-epochs 25` increases the max resolution from a Fourier radius of 12 pixels (`--l-start`) to 32 pix (`--l-end`) over the first 25 epochs of training.
    - Example training time course (1 V100 GPU)

        ```bash
        2022-01-20 18:00:59     # =====> Epoch: 1 Average gen loss = 0.8816, KLD = 0.8981, total loss = 0.8817; Finished in 1:13:20.758477
        2022-01-20 18:02:07     Using previous iteration poses
        2022-01-20 18:15:20     # =====> Epoch: 2 Average gen loss = 0.8825, KLD = 1.6127, total loss = 0.8826; Finished in 0:13:13.168348
        2022-01-20 18:16:27     Using previous iteration poses
        2022-01-20 18:29:41     # =====> Epoch: 3 Average gen loss = 0.8818, KLD = 1.8766, total loss = 0.8819; Finished in 0:13:14.378679
        2022-01-20 18:30:48     Using previous iteration poses
        2022-01-20 18:44:02     # =====> Epoch: 4 Average gen loss = 0.8811, KLD = 2.0323, total loss = 0.8812; Finished in 0:13:13.887047
        2022-01-20 18:45:09     Using previous iteration poses
        2022-01-20 18:58:23     # =====> Epoch: 5 Average gen loss = 0.8808, KLD = 2.1141, total loss = 0.8809; Finished in 0:13:14.298884
        2022-01-20 20:30:25     # =====> Epoch: 6 Average gen loss = 0.8783, KLD = 2.0354, total loss = 0.8784; Finished in 1:30:55.173547
        2022-01-20 20:31:46     Using previous iteration poses
        2022-01-20 20:45:00     # =====> Epoch: 7 Average gen loss = 0.878, KLD = 2.1416, total loss = 0.8782; Finished in 0:13:14.021982
        2022-01-20 20:46:20     Using previous iteration poses
        2022-01-20 20:59:35     # =====> Epoch: 8 Average gen loss = 0.8778, KLD = 2.1859, total loss = 0.8780; Finished in 0:13:14.906471
        2022-01-20 21:00:43     Using previous iteration poses
        2022-01-20 21:13:57     # =====> Epoch: 9 Average gen loss = 0.8776, KLD = 2.2222, total loss = 0.8778; Finished in 0:13:14.233966
        2022-01-20 21:15:04     Using previous iteration poses
        2022-01-20 21:28:19     # =====> Epoch: 10 Average gen loss = 0.8775, KLD = 2.2439, total loss = 0.8776; Finished in 0:13:14.907938
        2022-01-20 23:28:54     # =====> Epoch: 11 Average gen loss = 0.8769, KLD = 2.2428, total loss = 0.8771; Finished in 1:59:27.793799
        2022-01-20 23:30:01     Using previous iteration poses
        2022-01-20 23:43:15     # =====> Epoch: 12 Average gen loss = 0.8769, KLD = 2.3463, total loss = 0.8770; Finished in 0:13:14.289931
        2022-01-20 23:44:23     Using previous iteration poses
        2022-01-20 23:57:37     # =====> Epoch: 13 Average gen loss = 0.8767, KLD = 2.3692, total loss = 0.8769; Finished in 0:13:14.531232
        2022-01-20 23:58:45     Using previous iteration poses
        2022-01-21 00:11:59     # =====> Epoch: 14 Average gen loss = 0.8766, KLD = 2.3928, total loss = 0.8768; Finished in 0:13:14.821960
        2022-01-21 00:13:07     Using previous iteration poses
        2022-01-21 00:26:22     # =====> Epoch: 15 Average gen loss = 0.8765, KLD = 2.4063, total loss = 0.8767; Finished in 0:13:15.422771
        2022-01-21 02:58:58     # =====> Epoch: 16 Average gen loss = 0.8762, KLD = 2.3726, total loss = 0.8764; Finished in 2:31:28.195825
        2022-01-21 03:00:05     Using previous iteration poses
        2022-01-21 03:13:07     # =====> Epoch: 17 Average gen loss = 0.8762, KLD = 2.4672, total loss = 0.8764; Finished in 0:13:02.429271
        2022-01-21 03:14:14     Using previous iteration poses
        2022-01-21 03:27:17     # =====> Epoch: 18 Average gen loss = 0.876, KLD = 2.4911, total loss = 0.8762; Finished in 0:13:02.989886
        2022-01-21 03:28:24     Using previous iteration poses
        2022-01-21 03:41:27     # =====> Epoch: 19 Average gen loss = 0.876, KLD = 2.5077, total loss = 0.8762; Finished in 0:13:02.990354
        2022-01-21 03:42:34     Using previous iteration poses
        2022-01-21 03:55:37     # =====> Epoch: 20 Average gen loss = 0.8759, KLD = 2.5235, total loss = 0.8761; Finished in 0:13:02.651984
        2022-01-21 07:09:20     # =====> Epoch: 21 Average gen loss = 0.8756, KLD = 2.4737, total loss = 0.8758; Finished in 3:12:35.716970
        2022-01-21 07:10:27     Using previous iteration poses
        2022-01-21 07:23:30     # =====> Epoch: 22 Average gen loss = 0.8757, KLD = 2.5532, total loss = 0.8759; Finished in 0:13:02.696251
        2022-01-21 07:24:37     Using previous iteration poses
        2022-01-21 07:37:40     # =====> Epoch: 23 Average gen loss = 0.8756, KLD = 2.5753, total loss = 0.8758; Finished in 0:13:03.271111
        2022-01-21 07:38:47     Using previous iteration poses
        2022-01-21 07:51:52     # =====> Epoch: 24 Average gen loss = 0.8755, KLD = 2.5935, total loss = 0.8757; Finished in 0:13:05.296650
        2022-01-21 07:52:59     Using previous iteration poses
        2022-01-21 08:06:03     # =====> Epoch: 25 Average gen loss = 0.8755, KLD = 2.6057, total loss = 0.8757; Finished in 0:13:03.414230
        2022-01-21 12:13:51     # =====> Epoch: 26 Average gen loss = 0.8753, KLD = 2.5529, total loss = 0.8755; Finished in 4:06:41.267859
        2022-01-21 12:14:59     Using previous iteration poses
        2022-01-21 12:28:01     # =====> Epoch: 27 Average gen loss = 0.8753, KLD = 2.6321, total loss = 0.8755; Finished in 0:13:02.813657
        2022-01-21 12:29:09     Using previous iteration poses
        2022-01-21 12:42:11     # =====> Epoch: 28 Average gen loss = 0.8752, KLD = 2.6528, total loss = 0.8754; Finished in 0:13:02.511731
        2022-01-21 12:43:18     Using previous iteration poses
        2022-01-21 12:56:21     # =====> Epoch: 29 Average gen loss = 0.8752, KLD = 2.6635, total loss = 0.8754; Finished in 0:13:02.821824
        2022-01-21 12:57:28     Using previous iteration poses
        2022-01-21 13:10:31     # =====> Epoch: 30 Average gen loss = 0.8751, KLD = 2.6744, total loss = 0.8753; Finished in 0:13:02.737088
        ```

### Questions and contact

If you have any questions about the method or software, please file a GitHub issue:

[https://github.com/zhonge/cryodrgn/issues](https://github.com/zhonge/cryodrgn/issues)

Or post in the cryoDRGN Google Group: [https://groups.google.com/g/cryodrgn](https://groups.google.com/g/cryodrgn).


### Reference

CryoDRGN2 software was developed by Ellen Zhong & Adam Lerer with software support from Vineet Bansal.
If you find the ab initio tools in cryoDRGN useful, please cite:

Zhong, Lerer, Davis, Berger. ICCV 2021.

[https://openaccess.thecvf.com/content/ICCV2021/html/Zhong_CryoDRGN2_Ab_Initio_Neural_Reconstruction_of_3D_Protein_Structures_From_ICCV_2021_paper.html](https://openaccess.thecvf.com/content/ICCV2021/html/Zhong_CryoDRGN2_Ab_Initio_Neural_Reconstruction_of_3D_Protein_Structures_From_ICCV_2021_paper.html)
