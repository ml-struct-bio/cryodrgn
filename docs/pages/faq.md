# FAQ / Troubleshooting

- The expected heterogeneity is not present in the cryoDRGN reconstruction.
    - There may be too many junk particles or other imaging outliers in your dataset.
&nbsp;
- How do I update cryoDRGN versions?
    - See this section of the installation instructions
      [here](installation.md).

- How can I check on intermediate results during training?
    - See the run.log that is in the output directory. The `cryodrgn analyze` command may be run at any time on the
      intermediate epochs.

- Can I use the loss function or learning curve to tell when my training has converged?
    - The training curve is plotted in the Jupiter notebook. We usually use this to diagnose any training instabilities
      (spikes in the training curve). The loss function on the training set does not generally indicate convergence of
      the model.

- Do I need to use backproject_voxel? (Step 4 in the GitHub)
    - This step is used as a brief sanity check that the poses/CTF parameters have been parsed correctly.

- I am seeing weird artifacts (e.g. spherical ringing) in the reconstructed volumes.
    - This is most likely an issue with the input CTF parameters, or some type of aliasing from signal pre-processing,
      e.g. in signal-subtracted particle stacks. Please reach out if you run into these artifacts or file a github issue.

- Does cryoDRGN work with membrane protein complexes?
    - It does! However, we find that cryoDRGN can capture the heterogeneity of the micelle instead of the protein
      complex. You can try training on signal subtracted images.

- Does *particle crowdedness* affect the training result?
    - *Particle crowdedness* means signals from adjacent particles are included in cryoDRGN training.
      Particle crowding does affect the results (e.g. learning heterogeneity of neighboring particles).
      If you notice that this is a concern, you can reduce the reals-space windowing that is applied, by adding
      something like `--window .6` when using `cryodrgn train_vae`.

- Is cryoDRGN sensitive to orientation bias?
    - As long as sufficient orientational coverage exists, cryoDRGN should give reasonable results.

- How does one generally recognize a cluster as a junk cluster?  Are there any visible features in the UMAP that can
  be used?
    - This requires some digging into the results and images. You might be able to tell by looking at the maps or images
      in obvious cases. Sometimes you'll need to do follow-up analysis to validate, e.g. 2D classification.

- Do better separated clusters more likely mean distinct conformational states?  If a UMAP shows clusters that are
  close to one another (but with clear boundaries), can one interpret it as a set of conformational states that are very
  close to one another?
    - No, not in general. Sometimes there is "repetition" in the latent space where homogeneous states split into
      multiple clusters (maybe there is something subtle going on in the background). You may want to try out the new
      [landscape analysis tool](landscape_analysis.md)
      if you want to more quantitatively define conformational states.

More resources can be found by searching the [github issues](https://github.com/zhonge/cryodrgn/issues) page.
