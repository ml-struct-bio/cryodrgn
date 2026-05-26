"""Structured optional-arg groups for the dashboard command builder.

Mirrors ``add_argument_group`` titles and flags from:
``cryodrgn.commands.abinit``, ``train_vae``, ``train_nn``, ``train_dec``, and
``analyze*`` commands.
``--ctf`` and ``--datadir`` are rendered in the template (required row and
top-row dataset column), not listed here. Other required fields are template-only.

When CLI argument groups change, update this module and
``command_builder_cli_help`` together so the dashboard stays aligned.
"""

from __future__ import annotations

import copy
import os
from typing import Any

from cryodrgn.dashboard.command_builder_cli_help import (
    attach_help_to_groups,
    load_cli_help_maps,
)

# Widget kinds:
#   text, number, select — value optional; omit from command if empty
#   flag_true — append cli[0] when checked
#   flag_false — append cli[0] when unchecked (default-on feature disabled)
#   no_amp — checked = omit (use AMP); unchecked = --no-amp
#   norm2 — two text fields id + "_a" and id + "_b" (values after the flag)

Arg = dict[str, Any]

Group = dict[str, Any]

Schema = dict[str, list[Group]]


def default_outdir_for_command(cmd: str, workdir: str | None = None) -> str:
    """Default ``-o`` / ``--outdir``: ``001_<cmd>``, optionally under ``workdir``."""
    leaf = f"001_{cmd}"
    return os.path.join(workdir, leaf) if workdir else leaf


def arg_is_num_epochs(arg: object) -> bool:
    """True for ``-n`` / ``--num-epochs`` (total training epoch count)."""
    if not isinstance(arg, dict):
        return False
    cli = arg.get("cli")
    if not isinstance(cli, (list, tuple)):
        return False
    flags = {str(flag).lower() for flag in cli}
    return "-n" in flags or "--num-epochs" in flags


def arg_is_epoch_denominated(arg: object) -> bool:
    """True when any CLI flag name contains ``epoch`` (value is a count of epochs)."""
    if arg_is_num_epochs(arg):
        return False
    if not isinstance(arg, dict):
        return False
    cli = arg.get("cli")
    if not isinstance(cli, (list, tuple)):
        return False
    return any("epoch" in str(flag).lower() for flag in cli)


def arg_is_batch_size_denominated(arg: object) -> bool:
    """True when any CLI flag is a batch-size argument (value is a count of images)."""
    if not isinstance(arg, dict):
        return False
    cli = arg.get("cli")
    if not isinstance(cli, (list, tuple)):
        return False
    for flag in cli:
        s = str(flag).lower()
        if "batch-size" in s or s == "-b":
            return True
    return False


def _g(title: str, args: list[Arg], description: str = "") -> Group:
    return {"title": title, "args": args, "description": description}


# Short blurbs for GitHub Pages group cards
# (also used when ``description`` is omitted on a group).
COMMAND_BUILDER_GROUP_DESCRIPTIONS: dict[str, str] = {
    "Checkpoint & seed": "Restart from saved weights or poses; set RNG seed and verbosity.",
    "Dataset loading": "Particle subsets, datadir layout, lazy I/O, and RELION options.",
    "Logging": "Log frequency, checkpoint intervals, and timing diagnostics.",
    "Training parameters": "Epoch counts, batching, workers, and training schedule.",
    "Optimizers": "Learning rates, schedules, and optimizer hyperparameters.",
    "Masking": "Solvent mask generation, softness, and radial limits.",
    "Losses": "Reconstruction and regularization loss weights.",
    "Z / heterogeneity": "Latent dimensionality and continuous heterogeneity settings.",
    "Hypervolume": "Hypervolume decoder grid and related architecture options.",
    "Pretrain": "Warm-up and pretraining before full reconstruction.",
    "Pose search": "Ab initio pose search iterations and image subsampling.",
    "Normalization & analysis": "Input normalization and post-run analysis flags.",
    "I/O & logging": "Output layout, save cadence, and training logs.",
    "Tilt series parameters": "Tomography tilt geometry and per-tilt treatments.",
    "Pose SGD": "SGD steps for refining particle poses during training.",
    "Encoder Network": "VAE encoder depth, width, and activation choices.",
    "Decoder Network": "VAE decoder architecture and expressiveness.",
    "Network Architecture": "Layer sizes, activations, and backbone layout.",
    "Latent Variables": "Latent dimension, priors, and embedding behaviour.",
    "Backprojection parameters": "Batching, CTF weighting, half-maps, FSCs, and regularization.",
    "Tilt series options": "Tomography tilt geometry and dose per tilt.",
    "Run options": "Device, output paths, and skip flags for the analysis run.",
    "Volumes to generate": "Principal-component traversals, k-means samples, and per-axis density.",
    "Volume post-processing": "Apix, downsampling, filtering, cropping, and contrast for output MRCs.",
    "Extra arguments for volume generation": "Sketch size, Apix, downsampling, and volume-generation toggles.",
    "Extra arguments for mask generation": "Threshold, dilation, cosine edge, and custom mask for landscape volumes.",
    "Extra arguments for clustering": "Agglomerative linkage and cluster count for landscape grouping.",
    "Extra arguments for landscape visualization": "PCA and plot dimensionality for landscape figures.",
    "Volume generation arguments": "Training-volume count, downsampling, and skip flags for landscape_full.",
    "Volume mapping arguments": "Batch sizes, epochs, learning rate, and MLP architecture.",
    "Volume PC clustering arguments": "Neighbour count and Leiden resolution for volume PC clusters.",
}


def attach_group_descriptions(schema: Schema) -> None:
    """Fill ``description`` on each group for static command-builder cards."""
    for groups in schema.values():
        for g in groups:
            if not g.get("description"):
                g["description"] = COMMAND_BUILDER_GROUP_DESCRIPTIONS.get(
                    g["title"],
                    f"Optional flags for {g['title'].lower()}.",
                )


ABINIT_GROUPS: list[Group] = [
    _g(
        "Checkpoint & seed",
        [
            {"id": "ab_load", "cli": ["--load"], "w": "text"},
            {"id": "ab_load_poses", "cli": ["--load-poses"], "w": "text"},
            {"id": "ab_seed", "cli": ["--seed"], "w": "number"},
            {"id": "ab_verbose", "cli": ["-v"], "w": "flag_true"},
        ],
    ),
    _g(
        "Dataset loading",
        [
            {"id": "ab_ind", "cli": ["--ind"], "w": "text"},
            {"id": "ab_relion31", "cli": ["--relion31"], "w": "flag_true"},
            {"id": "ab_uninvert", "cli": ["--uninvert-data"], "w": "flag_true"},
            {"id": "ab_lazy", "cli": ["--lazy"], "w": "flag_true"},
            {
                "id": "ab_max_threads",
                "cli": ["--max-threads"],
                "w": "number",
                "placeholder": "16",
            },
        ],
    ),
    _g(
        "Logging",
        [
            {
                "id": "ab_log_interval",
                "cli": ["--log-interval"],
                "w": "number",
                "placeholder": "10000",
            },
            {
                "id": "ab_checkpoint",
                "cli": ["--checkpoint"],
                "w": "number",
                "placeholder": "5",
            },
            {"id": "ab_verbose_time", "cli": ["--verbose-time"], "w": "flag_true"},
        ],
    ),
    _g(
        "Training parameters",
        [
            {
                "id": "ab_n",
                "cli": ["-n", "--num-epochs"],
                "w": "number",
                "placeholder": "30",
            },
            {
                "id": "ab_epochs_pose_search",
                "cli": ["--epochs-pose-search"],
                "w": "number",
            },
            {
                "id": "ab_n_imgs_pose_search",
                "cli": ["--n-imgs-pose-search"],
                "w": "number",
            },
            {"id": "ab_epochs_sgd", "cli": ["--epochs-sgd"], "w": "number"},
            {
                "id": "ab_pose_only_phase",
                "cli": ["--pose-only-phase"],
                "w": "number",
                "placeholder": "0",
            },
            {"id": "ab_shuffle", "cli": ["--no-shuffle"], "w": "flag_false"},
            {
                "id": "ab_num_workers",
                "cli": ["--num-workers"],
                "w": "number",
                "placeholder": "2",
            },
            {
                "id": "ab_shuffler_size",
                "cli": ["--shuffler-size"],
                "w": "number",
                "placeholder": "32768",
            },
            {"id": "ab_multigpu", "cli": ["--multigpu"], "w": "flag_true"},
            {"id": "ab_use_amp", "w": "no_amp"},
            {
                "id": "ab_batch_size_hps",
                "cli": ["--batch-size-hps"],
                "w": "number",
                "placeholder": "16",
            },
            {
                "id": "ab_batch_size_known_poses",
                "cli": ["--batch-size-known-poses"],
                "w": "number",
                "placeholder": "64",
            },
            {
                "id": "ab_batch_size_sgd",
                "cli": ["--batch-size-sgd"],
                "w": "number",
                "placeholder": "128",
            },
        ],
    ),
    _g(
        "Optimizers",
        [
            {"id": "ab_lr", "cli": ["--lr"], "w": "text", "placeholder": "1e-4"},
            {
                "id": "ab_lr_pose_table",
                "cli": ["--lr-pose-table"],
                "w": "text",
                "placeholder": "1e-3",
            },
            {
                "id": "ab_lr_conf_table",
                "cli": ["--lr-conf-table"],
                "w": "text",
                "placeholder": "1e-2",
            },
            {
                "id": "ab_lr_conf_encoder",
                "cli": ["--lr-conf-encoder"],
                "w": "text",
                "placeholder": "1e-4",
            },
            {"id": "ab_wd", "cli": ["--wd"], "w": "text", "placeholder": "0"},
            {
                "id": "ab_hypervolume_optimizer_type",
                "cli": ["--hypervolume-optimizer-type"],
                "w": "select",
                "choices": ["adam"],
            },
            {
                "id": "ab_pose_table_optimizer_type",
                "cli": ["--pose-table-optimizer-type"],
                "w": "select",
                "choices": ["adam", "lbfgs"],
            },
            {
                "id": "ab_conf_table_optimizer_type",
                "cli": ["--conf-table-optimizer-type"],
                "w": "select",
                "choices": ["adam", "lbfgs"],
            },
            {
                "id": "ab_conf_encoder_optimizer_type",
                "cli": ["--conf-encoder-optimizer-type"],
                "w": "select",
                "choices": ["adam"],
            },
        ],
    ),
    _g(
        "Masking",
        [
            {
                "id": "ab_output_mask",
                "cli": ["--output-mask"],
                "w": "select",
                "choices": ["circ", "frequency_marching"],
            },
            {
                "id": "ab_add_one_frequency_every",
                "cli": ["--add-one-frequency-every"],
                "w": "number",
                "placeholder": "100000",
            },
            {
                "id": "ab_n_frequencies_per_epoch",
                "cli": ["--n-frequencies-per-epoch"],
                "w": "number",
                "placeholder": "10",
            },
            {"id": "ab_max_freq", "cli": ["--max-freq"], "w": "number"},
            {
                "id": "ab_window_radius_gt_real",
                "cli": ["--window-radius-gt-real"],
                "w": "text",
                "placeholder": "0.85",
            },
        ],
    ),
    _g(
        "Losses",
        [
            {
                "id": "ab_beta_conf",
                "cli": ["--beta-conf"],
                "w": "text",
                "placeholder": "0",
            },
            {
                "id": "ab_trans_l1_regularizer",
                "cli": ["--trans-l1-regularizer"],
                "w": "text",
                "placeholder": "0",
            },
            {
                "id": "ab_l2_smoothness_regularizer",
                "cli": ["--l2-smoothness-regularizer"],
                "w": "text",
                "placeholder": "0",
            },
        ],
    ),
    _g(
        "Z / heterogeneity",
        [
            {
                "id": "ab_variational_het",
                "cli": ["--variational-het"],
                "w": "flag_true",
            },
            {
                "id": "ab_std_z_init",
                "cli": ["--std-z-init"],
                "w": "text",
                "placeholder": "0.1",
            },
            {
                "id": "ab_use_conf_encoder",
                "cli": ["--use-conf-encoder"],
                "w": "flag_true",
            },
            {
                "id": "ab_depth_cnn",
                "cli": ["--depth-cnn"],
                "w": "number",
                "placeholder": "5",
            },
            {
                "id": "ab_channels_cnn",
                "cli": ["--channels-cnn"],
                "w": "number",
                "placeholder": "32",
            },
            {
                "id": "ab_kernel_size_cnn",
                "cli": ["--kernel-size-cnn"],
                "w": "number",
                "placeholder": "3",
            },
            {
                "id": "ab_resolution_encoder",
                "cli": ["--resolution-encoder"],
                "w": "number",
            },
        ],
    ),
    _g(
        "Hypervolume",
        [
            {
                "id": "ab_explicit_volume",
                "cli": ["--explicit-volume"],
                "w": "flag_true",
            },
            {"id": "ab_layers", "cli": ["--layers"], "w": "number", "placeholder": "3"},
            {"id": "ab_dim", "cli": ["--dim"], "w": "number", "placeholder": "256"},
            {
                "id": "ab_pe_type",
                "cli": ["--pe-type"],
                "w": "select",
                "choices": ["gaussian"],
            },
            {
                "id": "ab_pe_dim",
                "cli": ["--pe-dim"],
                "w": "number",
                "placeholder": "64",
            },
            {
                "id": "ab_feat_sigma",
                "cli": ["--feat-sigma"],
                "w": "text",
                "placeholder": "0.5",
            },
            {
                "id": "ab_hypervolume_domain",
                "cli": ["--hypervolume-domain"],
                "w": "select",
                "choices": ["hartley"],
            },
            {
                "id": "ab_pe_type_conf",
                "cli": ["--pe-type-conf"],
                "w": "select",
                "choices": ["", "geom"],
                "labels": ["(default None)", "geom"],
            },
            {"id": "ab_initial_conf", "cli": ["--initial-conf"], "w": "text"},
        ],
    ),
    _g(
        "Pretrain",
        [
            {
                "id": "ab_n_imgs_pretrain",
                "cli": ["--n-imgs-pretrain"],
                "w": "number",
                "placeholder": "10000",
            },
        ],
    ),
    _g(
        "Pose search",
        [
            {
                "id": "ab_l_start",
                "cli": ["--l-start"],
                "w": "number",
                "placeholder": "12",
            },
            {"id": "ab_l_end", "cli": ["--l-end"], "w": "number", "placeholder": "32"},
            {"id": "ab_niter", "cli": ["--niter"], "w": "number", "placeholder": "4"},
            {
                "id": "ab_t_extent",
                "cli": ["--t-extent"],
                "w": "text",
                "placeholder": "20.0",
            },
            {
                "id": "ab_t_ngrid",
                "cli": ["--t-ngrid"],
                "w": "number",
                "placeholder": "7",
            },
            {
                "id": "ab_t_xshift",
                "cli": ["--t-xshift"],
                "w": "text",
                "placeholder": "0.0",
            },
            {
                "id": "ab_t_yshift",
                "cli": ["--t-yshift"],
                "w": "text",
                "placeholder": "0.0",
            },
            {
                "id": "ab_no_trans_search_at_pose_search",
                "cli": ["--no-trans-search-at-pose-search"],
                "w": "flag_true",
            },
            {
                "id": "ab_nkeptposes",
                "cli": ["--nkeptposes"],
                "w": "number",
                "placeholder": "8",
            },
            {
                "id": "ab_base_healpy",
                "cli": ["--base-healpy"],
                "w": "number",
                "placeholder": "2",
            },
            {"id": "ab_no_trans", "cli": ["--no-trans"], "w": "flag_true"},
        ],
    ),
    _g(
        "Normalization & analysis",
        [
            {"id": "ab_norm", "cli": ["--norm"], "w": "norm2"},
            {"id": "ab_no_analysis", "cli": ["--no-analysis"], "w": "flag_false"},
        ],
    ),
]

TRAIN_VAE_GROUPS: list[Group] = [
    _g(
        "I/O & logging",
        [
            {"id": "vae_load", "cli": ["--load"], "w": "text"},
            {"id": "vae_no_analysis", "cli": ["--no-analysis"], "w": "flag_false"},
            {
                "id": "vae_checkpoint",
                "cli": ["--checkpoint"],
                "w": "number",
                "placeholder": "1",
            },
            {
                "id": "vae_log_interval",
                "cli": ["--log-interval"],
                "w": "number",
                "placeholder": "1000",
            },
            {"id": "vae_verbose", "cli": ["-v"], "w": "flag_true"},
            {"id": "vae_seed", "cli": ["--seed"], "w": "number"},
            {"id": "vae_shuffle_seed", "cli": ["--shuffle-seed"], "w": "number"},
        ],
    ),
    _g(
        "Dataset loading",
        [
            {"id": "vae_ind", "cli": ["--ind"], "w": "text"},
            {"id": "vae_uninvert", "cli": ["--uninvert-data"], "w": "flag_true"},
            {"id": "vae_window", "cli": ["--no-window"], "w": "flag_false"},
            {
                "id": "vae_window_r",
                "cli": ["--window-r"],
                "w": "text",
                "placeholder": "0.85",
            },
            {"id": "vae_lazy", "cli": ["--lazy"], "w": "flag_true"},
            {
                "id": "vae_shuffler_size",
                "cli": ["--shuffler-size"],
                "w": "number",
                "placeholder": "0",
            },
            {
                "id": "vae_num_workers",
                "cli": ["--num-workers"],
                "w": "number",
                "placeholder": "0",
            },
            {
                "id": "vae_max_threads",
                "cli": ["--max-threads"],
                "w": "number",
                "placeholder": "16",
            },
        ],
    ),
    _g(
        "Tilt series parameters",
        [
            {
                "id": "vae_ntilts",
                "cli": ["--ntilts"],
                "w": "number",
                "placeholder": "10",
            },
            {"id": "vae_random_tilts", "cli": ["--random-tilts"], "w": "flag_true"},
            {
                "id": "vae_t_emb_dim",
                "cli": ["--t-emb-dim"],
                "w": "number",
                "placeholder": "64",
            },
            {
                "id": "vae_tlayers",
                "cli": ["--tlayers"],
                "w": "number",
                "placeholder": "3",
            },
            {"id": "vae_tdim", "cli": ["--tdim"], "w": "number", "placeholder": "1024"},
            {"id": "vae_dose_per_tilt", "cli": ["-d", "--dose-per-tilt"], "w": "text"},
            {
                "id": "vae_angle_per_tilt",
                "cli": ["-a", "--angle-per-tilt"],
                "w": "text",
                "placeholder": "3",
            },
        ],
    ),
    _g(
        "Training parameters",
        [
            {
                "id": "vae_n",
                "cli": ["-n", "--num-epochs"],
                "w": "number",
                "placeholder": "20",
            },
            {
                "id": "vae_b",
                "cli": ["-b", "--batch-size"],
                "w": "number",
                "placeholder": "16",
            },
            {"id": "vae_wd", "cli": ["--wd"], "w": "text", "placeholder": "0"},
            {"id": "vae_lr", "cli": ["--lr"], "w": "text", "placeholder": "1e-4"},
            {"id": "vae_beta", "cli": ["--beta"], "w": "text"},
            {"id": "vae_beta_control", "cli": ["--beta-control"], "w": "text"},
            {"id": "vae_norm", "cli": ["--norm"], "w": "norm2"},
            {"id": "vae_use_amp", "w": "no_amp"},
            {"id": "vae_multigpu", "cli": ["--multigpu"], "w": "flag_true"},
        ],
    ),
    _g(
        "Pose SGD",
        [
            {"id": "vae_pose_sgd", "cli": ["--do-pose-sgd"], "w": "flag_true"},
            {
                "id": "vae_pretrain",
                "cli": ["--pretrain"],
                "w": "number",
                "placeholder": "1",
            },
            {
                "id": "vae_emb_type",
                "cli": ["--emb-type"],
                "w": "select",
                "choices": ["s2s2", "quat"],
            },
            {
                "id": "vae_pose_lr",
                "cli": ["--pose-lr"],
                "w": "text",
                "placeholder": "3e-4",
            },
        ],
    ),
    _g(
        "Encoder Network",
        [
            {
                "id": "vae_enc_layers",
                "cli": ["--enc-layers"],
                "w": "number",
                "placeholder": "3",
            },
            {
                "id": "vae_enc_dim",
                "cli": ["--enc-dim"],
                "w": "number",
                "placeholder": "1024",
            },
            {
                "id": "vae_encode_mode",
                "cli": ["--encode-mode"],
                "w": "select",
                "choices": ["conv", "resid", "mlp", "tilt"],
            },
            {"id": "vae_enc_mask", "cli": ["--enc-mask"], "w": "number"},
            {"id": "vae_use_real", "cli": ["--use-real"], "w": "flag_true"},
        ],
    ),
    _g(
        "Decoder Network",
        [
            {
                "id": "vae_dec_layers",
                "cli": ["--dec-layers"],
                "w": "number",
                "placeholder": "3",
            },
            {
                "id": "vae_dec_dim",
                "cli": ["--dec-dim"],
                "w": "number",
                "placeholder": "1024",
            },
            {
                "id": "vae_pe_type",
                "cli": ["--pe-type"],
                "w": "select",
                "choices": [
                    "geom_ft",
                    "geom_full",
                    "geom_lowf",
                    "geom_nohighf",
                    "linear_lowf",
                    "gaussian",
                    "none",
                ],
            },
            {
                "id": "vae_feat_sigma",
                "cli": ["--feat-sigma"],
                "w": "text",
                "placeholder": "0.5",
            },
            {"id": "vae_pe_dim", "cli": ["--pe-dim"], "w": "number"},
            {
                "id": "vae_domain",
                "cli": ["--domain"],
                "w": "select",
                "choices": ["hartley", "fourier"],
            },
            {
                "id": "vae_activation",
                "cli": ["--activation"],
                "w": "select",
                "choices": ["relu", "leaky_relu"],
            },
        ],
    ),
]

TRAIN_NN_GROUPS: list[Group] = [
    _g(
        "I/O & logging",
        [
            {"id": "nn_load", "cli": ["--load"], "w": "text"},
            {
                "id": "nn_checkpoint",
                "cli": ["--checkpoint"],
                "w": "number",
                "placeholder": "1",
            },
            {
                "id": "nn_log_interval",
                "cli": ["--log-interval"],
                "w": "number",
                "placeholder": "1000",
            },
            {"id": "nn_verbose", "cli": ["-v"], "w": "flag_true"},
            {"id": "nn_seed", "cli": ["--seed"], "w": "number"},
            {"id": "nn_shuffle_seed", "cli": ["--shuffle-seed"], "w": "number"},
        ],
    ),
    _g(
        "Dataset loading",
        [
            {"id": "nn_uninvert", "cli": ["--uninvert-data"], "w": "flag_true"},
            {"id": "nn_window", "cli": ["--no-window"], "w": "flag_false"},
            {
                "id": "nn_window_r",
                "cli": ["--window-r"],
                "w": "text",
                "placeholder": "0.85",
            },
            {"id": "nn_ind", "cli": ["--ind"], "w": "text"},
            {"id": "nn_lazy", "cli": ["--lazy"], "w": "flag_true"},
            {
                "id": "nn_shuffler_size",
                "cli": ["--shuffler-size"],
                "w": "number",
                "placeholder": "0",
            },
        ],
    ),
    _g(
        "Training parameters",
        [
            {
                "id": "nn_n",
                "cli": ["-n", "--num-epochs"],
                "w": "number",
                "placeholder": "20",
            },
            {
                "id": "nn_b",
                "cli": ["-b", "--batch-size"],
                "w": "number",
                "placeholder": "16",
            },
            {"id": "nn_wd", "cli": ["--wd"], "w": "text", "placeholder": "0"},
            {"id": "nn_lr", "cli": ["--lr"], "w": "text", "placeholder": "1e-4"},
            {"id": "nn_norm", "cli": ["--norm"], "w": "norm2"},
            {"id": "nn_use_amp", "w": "no_amp"},
            {"id": "nn_multigpu", "cli": ["--multigpu"], "w": "flag_true"},
        ],
    ),
    _g(
        "Pose SGD",
        [
            {"id": "nn_pose_sgd", "cli": ["--do-pose-sgd"], "w": "flag_true"},
            {
                "id": "nn_pretrain",
                "cli": ["--pretrain"],
                "w": "number",
                "placeholder": "5",
            },
            {
                "id": "nn_emb_type",
                "cli": ["--emb-type"],
                "w": "select",
                "choices": ["s2s2", "quat"],
            },
            {
                "id": "nn_pose_lr",
                "cli": ["--pose-lr"],
                "w": "text",
                "placeholder": "1e-4",
            },
        ],
    ),
    _g(
        "Network Architecture",
        [
            {"id": "nn_layers", "cli": ["--layers"], "w": "number", "placeholder": "3"},
            {"id": "nn_dim", "cli": ["--dim"], "w": "number", "placeholder": "1024"},
            {
                "id": "nn_l_extent",
                "cli": ["--l-extent"],
                "w": "text",
                "placeholder": "0.5",
            },
            {
                "id": "nn_pe_type",
                "cli": ["--pe-type"],
                "w": "select",
                "choices": [
                    "geom_ft",
                    "geom_full",
                    "geom_lowf",
                    "geom_nohighf",
                    "linear_lowf",
                    "gaussian",
                    "none",
                ],
            },
            {"id": "nn_pe_dim", "cli": ["--pe-dim"], "w": "number"},
            {
                "id": "nn_domain",
                "cli": ["--domain"],
                "w": "select",
                "choices": ["hartley", "fourier"],
            },
            {
                "id": "nn_activation",
                "cli": ["--activation"],
                "w": "select",
                "choices": ["relu", "leaky_relu"],
            },
            {
                "id": "nn_feat_sigma",
                "cli": ["--feat-sigma"],
                "w": "text",
                "placeholder": "0.5",
            },
        ],
    ),
]

TRAIN_DEC_GROUPS: list[Group] = [
    _g(
        "I/O & logging",
        [
            {"id": "dec_load", "cli": ["--load"], "w": "text"},
            {"id": "dec_no_analysis", "cli": ["--no-analysis"], "w": "flag_false"},
            {
                "id": "dec_checkpoint",
                "cli": ["--checkpoint"],
                "w": "number",
                "placeholder": "1",
            },
            {
                "id": "dec_log_interval",
                "cli": ["--log-interval"],
                "w": "number",
                "placeholder": "1000",
            },
            {"id": "dec_verbose", "cli": ["-v"], "w": "flag_true"},
            {"id": "dec_seed", "cli": ["--seed"], "w": "number"},
            {"id": "dec_shuffle_seed", "cli": ["--shuffle-seed"], "w": "number"},
        ],
    ),
    _g(
        "Latent Variables",
        [
            {"id": "dec_load_z", "cli": ["--load-z"], "w": "text"},
            {"id": "dec_z_lr", "cli": ["--z-lr"], "w": "text", "placeholder": "1e-4"},
            {
                "id": "dec_pretrain_z",
                "cli": ["--pretrain-z"],
                "w": "number",
                "placeholder": "0",
            },
        ],
    ),
    _g(
        "Dataset loading",
        [
            {"id": "dec_uninvert", "cli": ["--uninvert-data"], "w": "flag_true"},
            {"id": "dec_window", "cli": ["--no-window"], "w": "flag_false"},
            {
                "id": "dec_window_r",
                "cli": ["--window-r"],
                "w": "text",
                "placeholder": "0.85",
            },
            {"id": "dec_ind", "cli": ["--ind"], "w": "text"},
            {"id": "dec_lazy", "cli": ["--lazy"], "w": "flag_true"},
            {
                "id": "dec_shuffler_size",
                "cli": ["--shuffler-size"],
                "w": "number",
                "placeholder": "0",
            },
        ],
    ),
    _g(
        "Training parameters",
        [
            {
                "id": "dec_n",
                "cli": ["-n", "--num-epochs"],
                "w": "number",
                "placeholder": "20",
            },
            {
                "id": "dec_b",
                "cli": ["-b", "--batch-size"],
                "w": "number",
                "placeholder": "8",
            },
            {"id": "dec_wd", "cli": ["--wd"], "w": "text", "placeholder": "0"},
            {"id": "dec_lr", "cli": ["--lr"], "w": "text", "placeholder": "1e-4"},
            {"id": "dec_norm", "cli": ["--norm"], "w": "norm2"},
            {"id": "dec_use_amp", "w": "no_amp"},
            {"id": "dec_multigpu", "cli": ["--multigpu"], "w": "flag_true"},
        ],
    ),
    _g(
        "Pose SGD",
        [
            {"id": "dec_pose_sgd", "cli": ["--do-pose-sgd"], "w": "flag_true"},
            {
                "id": "dec_pretrain_pose",
                "cli": ["--pretrain-pose"],
                "w": "number",
                "placeholder": "5",
            },
            {
                "id": "dec_emb_type",
                "cli": ["--emb-type"],
                "w": "select",
                "choices": ["s2s2", "quat"],
            },
            {
                "id": "dec_pose_lr",
                "cli": ["--pose-lr"],
                "w": "text",
                "placeholder": "1e-4",
            },
        ],
    ),
    _g(
        "Network Architecture",
        [
            {
                "id": "dec_layers",
                "cli": ["--layers"],
                "w": "number",
                "placeholder": "3",
            },
            {"id": "dec_dim", "cli": ["--dim"], "w": "number", "placeholder": "1024"},
            {
                "id": "dec_l_extent",
                "cli": ["--l-extent"],
                "w": "text",
                "placeholder": "0.5",
            },
            {
                "id": "dec_pe_type",
                "cli": ["--pe-type"],
                "w": "select",
                "choices": [
                    "geom_ft",
                    "geom_full",
                    "geom_lowf",
                    "geom_nohighf",
                    "linear_lowf",
                    "gaussian",
                    "none",
                ],
            },
            {"id": "dec_pe_dim", "cli": ["--pe-dim"], "w": "number"},
            {
                "id": "dec_domain",
                "cli": ["--domain"],
                "w": "select",
                "choices": ["hartley", "fourier"],
            },
            {
                "id": "dec_activation",
                "cli": ["--activation"],
                "w": "select",
                "choices": ["relu", "leaky_relu"],
            },
            {
                "id": "dec_feat_sigma",
                "cli": ["--feat-sigma"],
                "w": "text",
                "placeholder": "0.5",
            },
        ],
    ),
]

BACKPROJECT_VOXEL_GROUPS: list[Group] = [
    _g(
        "Dataset loading",
        [
            {"id": "bpv_uninvert", "cli": ["--uninvert-data"], "w": "flag_true"},
            {"id": "bpv_lazy", "cli": ["--lazy"], "w": "flag_true"},
            {"id": "bpv_ind", "cli": ["--ind"], "w": "text"},
            {"id": "bpv_first", "cli": ["--first"], "w": "number"},
        ],
    ),
    _g(
        "Backprojection parameters",
        [
            {
                "id": "bpv_half_maps",
                "cli": ["--no-half-maps"],
                "w": "flag_false",
            },
            {
                "id": "bpv_fsc_vals",
                "cli": ["--no-fsc-vals"],
                "w": "flag_false",
            },
            {
                "id": "bpv_batch_size",
                "cli": ["-b", "--batch-size"],
                "w": "number",
                "placeholder": "1000",
            },
            {
                "id": "bpv_ctf_alg",
                "cli": ["--ctf-alg"],
                "w": "select",
                "choices": ["flip", "mul"],
                "default": "mul",
                "help": "CTF algorithm: phase flip (flip) or multiply (mul).",
            },
            {
                "id": "bpv_reg_weight",
                "cli": ["--reg-weight"],
                "w": "text",
                "placeholder": "1.0",
            },
            {
                "id": "bpv_output_sumcount",
                "cli": ["--output-sumcount"],
                "w": "flag_true",
            },
            {
                "id": "bpv_log_interval",
                "cli": ["--log-interval"],
                "w": "text",
                "placeholder": "5000",
            },
        ],
    ),
    _g(
        "Tilt series parameters",
        [
            {"id": "bpv_tilt", "cli": ["--tilt"], "w": "flag_true"},
            {
                "id": "bpv_ntilts",
                "cli": ["--ntilts"],
                "w": "number",
                "placeholder": "10",
            },
            {"id": "bpv_force_ntilts", "cli": ["--force-ntilts"], "w": "flag_true"},
            {"id": "bpv_dose_per_tilt", "cli": ["-d", "--dose-per-tilt"], "w": "text"},
            {
                "id": "bpv_angle_per_tilt",
                "cli": ["-a", "--angle-per-tilt"],
                "w": "text",
                "placeholder": "3",
            },
        ],
    ),
]

# Primary-band group titles surfaced in the GitHub Pages dataset / run-parameter pair.
COMMAND_BUILDER_PRIMARY_RUN_GROUP_TITLES: frozenset[str] = frozenset(
    {
        "Training parameters",
        "Backprojection parameters",
        "Run options",
        "Volumes to generate",
        "Extra arguments for volume generation",
        "Volume generation arguments",
        "Volume mapping arguments",
    }
)


def _remap_group_ids(
    groups: list[Group], from_prefix: str, to_prefix: str
) -> list[Group]:
    out: list[Group] = []
    for g in groups:
        ng = copy.deepcopy(g)
        for a in ng["args"]:
            aid = a.get("id", "")
            if aid.startswith(from_prefix):
                a["id"] = to_prefix + aid[len(from_prefix) :]
        out.append(ng)
    return out


def _arg_in_help_map(a: Arg, help_map: dict[str, str]) -> bool:
    w = a.get("w")
    if w == "no_amp":
        return "--no-amp" in help_map
    cli = a.get("cli") or []
    return bool(cli) and any(c in help_map for c in cli)


def _filter_groups_to_help_map(
    groups: list[Group], help_map: dict[str, str]
) -> list[Group]:
    filtered: list[Group] = []
    for g in groups:
        args = [a for a in g["args"] if _arg_in_help_map(a, help_map)]
        if args:
            filtered.append({**g, "args": args})
    return filtered


def _build_abinit_het_old_groups() -> list[Group]:
    """Schema for deprecated ``abinit_het_old`` (VAE + hierarchical pose search)."""
    groups = _remap_group_ids(TRAIN_VAE_GROUPS, "vae_", "ahet_")
    groups = [g for g in groups if g["title"] != "Pose SGD"]
    pose_search = _g(
        "Pose Search parameters",
        [
            {
                "id": "ahet_l_start",
                "cli": ["--l-start"],
                "w": "number",
                "placeholder": "12",
            },
            {
                "id": "ahet_l_end",
                "cli": ["--l-end"],
                "w": "number",
                "placeholder": "32",
            },
            {"id": "ahet_niter", "cli": ["--niter"], "w": "number", "placeholder": "4"},
            {
                "id": "ahet_t_extent",
                "cli": ["--t-extent"],
                "w": "text",
                "placeholder": "10",
            },
            {
                "id": "ahet_t_ngrid",
                "cli": ["--t-ngrid"],
                "w": "number",
                "placeholder": "7",
            },
            {
                "id": "ahet_t_xshift",
                "cli": ["--t-xshift"],
                "w": "text",
                "placeholder": "0",
            },
            {
                "id": "ahet_t_yshift",
                "cli": ["--t-yshift"],
                "w": "text",
                "placeholder": "0",
            },
            {
                "id": "ahet_pretrain_ps",
                "cli": ["--pretrain"],
                "w": "number",
                "placeholder": "10000",
            },
            {
                "id": "ahet_ps_freq",
                "cli": ["--ps-freq"],
                "w": "number",
                "placeholder": "5",
            },
            {
                "id": "ahet_nkeptposes",
                "cli": ["--nkeptposes"],
                "w": "number",
                "placeholder": "8",
            },
            {
                "id": "ahet_base_healpy",
                "cli": ["--base-healpy"],
                "w": "number",
                "placeholder": "2",
            },
            {
                "id": "ahet_pose_model_update_freq",
                "cli": ["--pose-model-update-freq"],
                "w": "number",
            },
        ],
    )
    for i, g in enumerate(groups):
        if g["title"] == "Tilt series parameters":
            g["title"] = "Tilt series"
            g["args"].append(
                {"id": "ahet_enc_only", "cli": ["--enc-only"], "w": "flag_true"}
            )
        elif g["title"] == "I/O & logging":
            g["args"].insert(
                1, {"id": "ahet_load_poses", "cli": ["--load-poses"], "w": "text"}
            )
        elif g["title"] == "Training parameters":
            g["args"].extend(
                [
                    {
                        "id": "ahet_equivariance",
                        "cli": ["--equivariance"],
                        "w": "text",
                    },
                    {
                        "id": "ahet_eq_start_it",
                        "cli": ["--eq-start-it"],
                        "w": "number",
                    },
                    {"id": "ahet_eq_end_it", "cli": ["--eq-end-it"], "w": "number"},
                    {
                        "id": "ahet_l_ramp_epochs",
                        "cli": ["--l-ramp-epochs"],
                        "w": "number",
                        "placeholder": "0",
                    },
                    {
                        "id": "ahet_l_ramp_model",
                        "cli": ["--l-ramp-model"],
                        "w": "number",
                        "placeholder": "0",
                    },
                    {
                        "id": "ahet_reset_model_every",
                        "cli": ["--reset-model-every"],
                        "w": "number",
                    },
                    {
                        "id": "ahet_reset_optim_every",
                        "cli": ["--reset-optim-every"],
                        "w": "number",
                    },
                    {
                        "id": "ahet_reset_optim_after_pretrain",
                        "cli": ["--reset-optim-after-pretrain"],
                        "w": "number",
                    },
                ]
            )
            groups.insert(i + 1, pose_search)
            break
    hm = load_cli_help_maps().get("abinit_het_old", {})
    return _filter_groups_to_help_map(groups, hm)


def _build_abinit_homo_old_groups() -> list[Group]:
    """Schema for deprecated ``cryodrgn abinit_homo_old`` (homo NN + pose search)."""
    groups = _remap_group_ids(TRAIN_NN_GROUPS, "nn_", "ahom_")
    groups = [g for g in groups if g["title"] != "Pose SGD"]
    tilt = _g(
        "Tilt series",
        [
            {"id": "ahom_tilt", "cli": ["--tilt"], "w": "text"},
            {
                "id": "ahom_tilt_deg",
                "cli": ["--tilt-deg"],
                "w": "text",
                "placeholder": "45",
            },
        ],
    )
    pose_search = _g(
        "Pose search parameters",
        [
            {
                "id": "ahom_l_start",
                "cli": ["--l-start"],
                "w": "number",
                "placeholder": "12",
            },
            {
                "id": "ahom_l_end",
                "cli": ["--l-end"],
                "w": "number",
                "placeholder": "32",
            },
            {"id": "ahom_niter", "cli": ["--niter"], "w": "number", "placeholder": "4"},
            {
                "id": "ahom_l_ramp_epochs",
                "cli": ["--l-ramp-epochs"],
                "w": "number",
                "placeholder": "25",
            },
            {"id": "ahom_probabilistic", "cli": ["--probabilistic"], "w": "flag_true"},
            {
                "id": "ahom_nkeptposes",
                "cli": ["--nkeptposes"],
                "w": "number",
                "placeholder": "8",
            },
            {
                "id": "ahom_base_healpy",
                "cli": ["--base-healpy"],
                "w": "number",
                "placeholder": "2",
            },
            {
                "id": "ahom_pose_model_update_freq",
                "cli": ["--pose-model-update-freq"],
                "w": "number",
            },
        ],
    )
    homo_train_extra: list[Arg] = [
        {
            "id": "ahom_t_extent",
            "cli": ["--t-extent"],
            "w": "text",
            "placeholder": "10",
        },
        {"id": "ahom_t_ngrid", "cli": ["--t-ngrid"], "w": "number", "placeholder": "7"},
        {"id": "ahom_t_xshift", "cli": ["--t-xshift"], "w": "text", "placeholder": "0"},
        {"id": "ahom_t_yshift", "cli": ["--t-yshift"], "w": "text", "placeholder": "0"},
        {"id": "ahom_no_trans", "cli": ["--no-trans"], "w": "flag_true"},
        {
            "id": "ahom_pretrain_train",
            "cli": ["--pretrain"],
            "w": "number",
            "placeholder": "10000",
        },
        {"id": "ahom_ps_freq", "cli": ["--ps-freq"], "w": "number", "placeholder": "5"},
    ]
    for i, g in enumerate(groups):
        if g["title"] == "I/O & logging":
            g["args"].insert(
                1, {"id": "ahom_load_poses", "cli": ["--load-poses"], "w": "text"}
            )
        elif g["title"] == "Dataset loading":
            groups.insert(i + 1, tilt)
        elif g["title"] == "Training parameters":
            g["args"] = homo_train_extra + g["args"]
            groups.insert(i + 1, pose_search)
            break
    hm = load_cli_help_maps().get("abinit_homo_old", {})
    return _filter_groups_to_help_map(groups, hm)


ABINIT_HET_OLD_GROUPS: list[Group] = _build_abinit_het_old_groups()
ABINIT_HOMO_OLD_GROUPS: list[Group] = _build_abinit_homo_old_groups()

ANALYZE_GROUPS: list[Group] = [
    _g(
        "Run options",
        [
            {"id": "ana_device", "cli": ["--device"], "w": "number"},
            {"id": "ana_skip_vol", "cli": ["--skip-vol"], "w": "flag_true"},
            {"id": "ana_skip_umap", "cli": ["--skip-umap"], "w": "flag_true"},
        ],
    ),
    _g(
        "Volumes to generate",
        [
            {"id": "ana_pc", "cli": ["--pc"], "w": "number", "placeholder": "2"},
            {
                "id": "ana_n_per_pc",
                "cli": ["--n-per-pc"],
                "w": "number",
                "placeholder": "10",
            },
            {
                "id": "ana_ksample",
                "cli": ["--ksample"],
                "w": "number",
                "placeholder": "20",
            },
        ],
    ),
    _g(
        "Volume post-processing",
        [
            {"id": "ana_apix", "cli": ["--Apix"], "w": "text"},
            {"id": "ana_flip", "cli": ["--flip"], "w": "flag_true"},
            {"id": "ana_invert", "cli": ["--invert"], "w": "flag_true"},
            {
                "id": "ana_downsample",
                "cli": ["-d", "--downsample"],
                "w": "number",
            },
            {"id": "ana_low_pass", "cli": ["--low-pass"], "w": "text"},
            {"id": "ana_crop", "cli": ["--crop"], "w": "number"},
            {
                "id": "ana_vol_start_index",
                "cli": ["--vol-start-index"],
                "w": "number",
                "placeholder": "1",
            },
        ],
    ),
]

ANALYZE_LANDSCAPE_GROUPS: list[Group] = [
    _g(
        "Run options",
        [
            {"id": "alsc_device", "cli": ["--device"], "w": "number"},
            {"id": "alsc_multigpu", "cli": ["--multigpu"], "w": "flag_true"},
            {"id": "alsc_skip_umap", "cli": ["--skip-umap"], "w": "flag_true"},
            {"id": "alsc_vol_ind", "cli": ["--vol-ind"], "w": "text"},
        ],
    ),
    _g(
        "Extra arguments for volume generation",
        [
            {
                "id": "alsc_sketch_size",
                "cli": ["-N", "--sketch-size"],
                "w": "number",
                "placeholder": "1000",
            },
            {
                "id": "alsc_apix",
                "cli": ["--Apix"],
                "w": "text",
                "placeholder": "1",
            },
            {"id": "alsc_flip", "cli": ["--flip"], "w": "flag_true"},
            {
                "id": "alsc_downsample",
                "cli": ["-d", "--downsample"],
                "w": "number",
                "placeholder": "128",
            },
            {"id": "alsc_skip_vol", "cli": ["--skip-vol"], "w": "flag_true"},
            {
                "id": "alsc_vol_start_index",
                "cli": ["--vol-start-index"],
                "w": "number",
                "placeholder": "1",
            },
        ],
    ),
    _g(
        "Extra arguments for mask generation",
        [
            {"id": "alsc_thresh", "cli": ["--thresh"], "w": "text"},
            {
                "id": "alsc_dilate",
                "cli": ["--dilate"],
                "w": "number",
                "placeholder": "5",
            },
            {
                "id": "alsc_cosine_edge",
                "cli": ["--cosine-edge"],
                "w": "number",
                "placeholder": "0",
            },
            {"id": "alsc_mask", "cli": ["--mask"], "w": "text"},
        ],
    ),
    _g(
        "Extra arguments for clustering",
        [
            {
                "id": "alsc_linkage",
                "cli": ["--linkage"],
                "w": "text",
                "placeholder": "average",
            },
            {
                "id": "alsc_n_clusters",
                "cli": ["-M"],
                "w": "number",
                "placeholder": "10",
            },
        ],
    ),
    _g(
        "Extra arguments for landscape visualization",
        [
            {
                "id": "alsc_pc_dim",
                "cli": ["--pc-dim"],
                "w": "number",
                "placeholder": "20",
            },
            {
                "id": "alsc_plot_dim",
                "cli": ["--plot-dim"],
                "w": "number",
                "placeholder": "5",
            },
        ],
    ),
]

ANALYZE_LANDSCAPE_FULL_GROUPS: list[Group] = [
    _g(
        "Run options",
        [
            {"id": "alfull_device", "cli": ["--device"], "w": "number"},
            {
                "id": "alfull_landscape_dir",
                "cli": ["--landscape-dir"],
                "w": "text",
            },
            {
                "id": "alfull_seed",
                "cli": ["--seed"],
                "w": "number",
                "placeholder": "0",
            },
        ],
    ),
    _g(
        "Volume generation arguments",
        [
            {
                "id": "alfull_training_volumes",
                "cli": ["-N", "--training-volumes"],
                "w": "number",
                "placeholder": "10000",
            },
            {"id": "alfull_flip", "cli": ["--flip"], "w": "flag_true"},
            {
                "id": "alfull_downsample",
                "cli": ["-d", "--downsample"],
                "w": "number",
                "placeholder": "128",
            },
            {"id": "alfull_skip_vol", "cli": ["--skip-vol"], "w": "flag_true"},
        ],
    ),
    _g(
        "Volume mapping arguments",
        [
            {
                "id": "alfull_batch_size",
                "cli": ["--batch-size"],
                "w": "number",
                "placeholder": "64",
            },
            {
                "id": "alfull_test_batch_size",
                "cli": ["--test-batch-size"],
                "w": "number",
                "placeholder": "1000",
            },
            {
                "id": "alfull_epochs",
                "cli": ["--epochs"],
                "w": "number",
                "placeholder": "200",
            },
            {
                "id": "alfull_lr",
                "cli": ["--lr"],
                "w": "text",
                "placeholder": "1e-4",
            },
            {
                "id": "alfull_dim",
                "cli": ["--dim"],
                "w": "number",
                "placeholder": "512",
            },
            {
                "id": "alfull_layers",
                "cli": ["--layers"],
                "w": "number",
                "placeholder": "3",
            },
        ],
    ),
    _g(
        "Volume PC clustering arguments",
        [
            {
                "id": "alfull_num_neighbors",
                "cli": ["--num-neighbors"],
                "w": "number",
                "placeholder": "50",
            },
            {
                "id": "alfull_resolution",
                "cli": ["--resolution"],
                "w": "text",
                "placeholder": "1.5",
            },
        ],
    ),
]

COMMAND_BUILDER_COMMAND_KEYS: tuple[str, ...] = (
    "abinit",
    "abinit_het_old",
    "abinit_homo_old",
    "train_vae",
    "train_nn",
    "train_dec",
    "backproject_voxel",
    "analyze",
    "analyze_landscape",
    "analyze_landscape_full",
)

COMMAND_BUILDER_SCHEMA: Schema = {
    "abinit": ABINIT_GROUPS,
    "abinit_het_old": ABINIT_HET_OLD_GROUPS,
    "abinit_homo_old": ABINIT_HOMO_OLD_GROUPS,
    "train_vae": TRAIN_VAE_GROUPS,
    "train_nn": TRAIN_NN_GROUPS,
    "train_dec": TRAIN_DEC_GROUPS,
    "backproject_voxel": BACKPROJECT_VOXEL_GROUPS,
    "analyze": ANALYZE_GROUPS,
    "analyze_landscape": ANALYZE_LANDSCAPE_GROUPS,
    "analyze_landscape_full": ANALYZE_LANDSCAPE_FULL_GROUPS,
}

_cli_help = load_cli_help_maps()
attach_help_to_groups(_cli_help.get("abinit", {}), ABINIT_GROUPS)
attach_help_to_groups(_cli_help.get("abinit_het_old", {}), ABINIT_HET_OLD_GROUPS)
attach_help_to_groups(_cli_help.get("abinit_homo_old", {}), ABINIT_HOMO_OLD_GROUPS)
attach_help_to_groups(_cli_help.get("train_vae", {}), TRAIN_VAE_GROUPS)
attach_help_to_groups(_cli_help.get("train_nn", {}), TRAIN_NN_GROUPS)
attach_help_to_groups(_cli_help.get("train_dec", {}), TRAIN_DEC_GROUPS)
attach_help_to_groups(_cli_help.get("backproject_voxel", {}), BACKPROJECT_VOXEL_GROUPS)
attach_help_to_groups(_cli_help.get("analyze", {}), ANALYZE_GROUPS)
attach_help_to_groups(_cli_help.get("analyze_landscape", {}), ANALYZE_LANDSCAPE_GROUPS)
attach_help_to_groups(
    _cli_help.get("analyze_landscape_full", {}),
    ANALYZE_LANDSCAPE_FULL_GROUPS,
)
attach_group_descriptions(COMMAND_BUILDER_SCHEMA)


def _required_field_titles(
    hm: dict[str, str],
    pairs: dict[str, str | tuple[str, ...]],
) -> dict[str, str]:
    out: dict[str, str] = {}
    for elid, keys in pairs.items():
        seq: tuple[str, ...] = (keys,) if isinstance(keys, str) else keys
        for k in seq:
            t = hm.get(k)
            if t:
                out[elid] = t
                break
    return out


def _build_required_field_titles() -> dict[str, str]:
    r: dict[str, str] = {}
    r.update(
        _required_field_titles(
            _cli_help.get("abinit", {}),
            {
                "ab_particles": "particles",
                "ab_out": ("-o", "--outdir"),
                "ab_zdim": "--zdim",
            },
        ),
    )
    r.update(
        _required_field_titles(
            _cli_help.get("abinit_het_old", {}),
            {
                "ahet_particles": "particles",
                "ahet_out": ("-o", "--outdir"),
                "ahet_zdim": "--zdim",
            },
        ),
    )
    r.update(
        _required_field_titles(
            _cli_help.get("abinit_homo_old", {}),
            {
                "ahom_particles": "particles",
                "ahom_out": ("-o", "--outdir"),
            },
        ),
    )
    r.update(
        _required_field_titles(
            _cli_help.get("train_vae", {}),
            {
                "vae_particles": "particles",
                "vae_out": ("-o", "--outdir"),
                "vae_poses": "--poses",
                "vae_zdim": "--zdim",
            },
        ),
    )
    r.update(
        _required_field_titles(
            _cli_help.get("train_nn", {}),
            {
                "nn_particles": "particles",
                "nn_out": ("-o", "--outdir"),
                "nn_poses": "--poses",
            },
        ),
    )
    r.update(
        _required_field_titles(
            _cli_help.get("train_dec", {}),
            {
                "dec_particles": "particles",
                "dec_out": ("-o", "--outdir"),
                "dec_poses": "--poses",
                "dec_zdim": "--zdim",
            },
        ),
    )
    r.update(
        _required_field_titles(
            _cli_help.get("backproject_voxel", {}),
            {
                "bpv_particles": "particles",
                "bpv_out": ("-o", "--outdir"),
                "bpv_poses": "--poses",
            },
        ),
    )
    r.update(
        _required_field_titles(
            _cli_help.get("analyze", {}),
            {
                "ana_workdir": "workdir",
                "ana_epoch": "epoch",
                "ana_out": ("-o", "--outdir"),
            },
        ),
    )
    r.update(
        _required_field_titles(
            _cli_help.get("analyze_landscape", {}),
            {
                "alsc_workdir": "workdir",
                "alsc_epoch": "epoch",
                "alsc_out": ("-o", "--outdir"),
            },
        ),
    )
    r.update(
        _required_field_titles(
            _cli_help.get("analyze_landscape_full", {}),
            {
                "alfull_workdir": "workdir",
                "alfull_epoch": "epoch",
                "alfull_out": ("-o", "--outdir"),
                "alfull_landscape_dir": "--landscape-dir",
            },
        ),
    )
    for prefix, cmd in (
        ("ab", "abinit"),
        ("ahet", "abinit_het_old"),
        ("ahom", "abinit_homo_old"),
        ("vae", "train_vae"),
        ("nn", "train_nn"),
        ("dec", "train_dec"),
        ("bpv", "backproject_voxel"),
    ):
        hm = _cli_help.get(cmd, {})
        if "--ctf" in hm:
            r[f"{prefix}_ctf"] = hm["--ctf"]
        if "--datadir" in hm:
            r[f"{prefix}_datadir"] = hm["--datadir"]
    return r


COMMAND_BUILDER_REQUIRED_FIELD_TITLES: dict[str, str] = _build_required_field_titles()
