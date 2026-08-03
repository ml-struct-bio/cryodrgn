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
    load_cli_defaults_maps,
    load_cli_help_maps,
    resolved_help_for_flag,
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


_EPOCH_COUNT_FLAGS = frozenset({"--pose-only-phase"})


def arg_is_epoch_denominated(arg: object) -> bool:
    """True when any CLI flag name contains ``epoch`` (value is a count of epochs)."""
    if arg_is_num_epochs(arg):
        return False
    if not isinstance(arg, dict):
        return False
    cli = arg.get("cli")
    if not isinstance(cli, (list, tuple)):
        return False
    for flag in cli:
        s = str(flag).lower()
        if s in _EPOCH_COUNT_FLAGS or "epoch" in s:
            return True
    return False


def arg_has_unit_suffix(arg: object) -> bool:
    """True when the form shows an ``epochs`` / ``images`` unit beside the value."""
    return (
        arg_is_num_epochs(arg)
        or arg_is_epoch_denominated(arg)
        or arg_is_batch_size_denominated(arg)
    )


_SUPPRESS_DISPLAY_NAME_FLAGS = frozenset({"--max-threads"})


def arg_show_display_name(arg: object) -> bool:
    """Whether to show the short human label beside the CLI flag."""
    if not isinstance(arg, dict):
        return False
    cli = arg.get("cli")
    if isinstance(cli, (list, tuple)):
        flags = {str(flag).lower() for flag in cli}
        if flags & _SUPPRESS_DISPLAY_NAME_FLAGS:
            return False
        if "-b" in flags:
            return True
    return not arg_has_unit_suffix(arg)


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


def batch_size_arg_ids_for_command(cmd: str) -> list[str]:
    """Return form field ids for batch-size CLI args on ``cmd``."""
    ids: list[str] = []
    for g in COMMAND_BUILDER_SCHEMA.get(cmd, []):
        for a in g.get("args", []):
            if a.get("w") == "number" and arg_is_batch_size_denominated(a):
                ids.append(str(a["id"]))
    return ids


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


# Argument factories for common patterns
def _text(id: str, cli: list[str]) -> Arg:
    """Text widget argument."""
    return {"id": id, "cli": cli, "w": "text"}


def _num(id: str, cli: list[str], placeholder: str | None = None) -> Arg:
    """Number widget argument with optional placeholder."""
    arg: Arg = {"id": id, "cli": cli, "w": "number"}
    if placeholder is not None:
        arg["placeholder"] = placeholder
    return arg


def _flag(id: str, cli: list[str]) -> Arg:
    """Flag_true widget argument."""
    return {"id": id, "cli": cli, "w": "flag_true"}


def _flag_false(id: str, cli: list[str]) -> Arg:
    """Flag_false widget argument."""
    return {"id": id, "cli": cli, "w": "flag_false"}


def _no_amp(id: str) -> Arg:
    """no_amp widget argument (checked = omit, unchecked = --no-amp)."""
    return {"id": id, "w": "no_amp"}


def _select(id: str, cli: list[str], choices: list[str], **extra: Any) -> Arg:
    """Select widget argument with choices and optional extra fields."""
    arg: Arg = {"id": id, "cli": cli, "w": "select", "choices": choices}
    arg.update(extra)
    return arg


# Common argument groups shared across commands
def _io_logging_group(prefix: str) -> Group:
    """I/O & logging group shared across training commands."""
    return _g(
        "I/O & logging",
        [
            _text(f"{prefix}_load", ["--load"]),
            _flag_false(f"{prefix}_no_analysis", ["--no-analysis"]),
            _num(f"{prefix}_checkpoint", ["--checkpoint"], "1"),
            _num(f"{prefix}_log_interval", ["--log-interval"], "1000"),
            _flag(f"{prefix}_verbose", ["-v"]),
            _num(f"{prefix}_seed", ["--seed"]),
            _num(f"{prefix}_shuffle_seed", ["--shuffle-seed"]),
        ],
    )


def _dataset_loading_group(prefix: str) -> Group:
    """Dataset loading group shared across training commands."""
    return _g(
        "Dataset loading",
        [
            _flag(f"{prefix}_uninvert", ["--uninvert-data"]),
            _flag_false(f"{prefix}_window", ["--no-window"]),
            _num(f"{prefix}_window_r", ["--window-r"], "0.85"),
            _text(f"{prefix}_ind", ["--ind"]),
            _flag(f"{prefix}_lazy", ["--lazy"]),
            _num(f"{prefix}_shuffler_size", ["--shuffler-size"], "0"),
        ],
    )


def _pose_sgd_group(prefix: str) -> Group:
    """Pose SGD group shared across training commands."""
    return _g(
        "Pose SGD",
        [
            _flag(f"{prefix}_pose_sgd", ["--do-pose-sgd"]),
            _num(f"{prefix}_pretrain", ["--pretrain"], "1"),
            _select(f"{prefix}_emb_type", ["--emb-type"], ["s2s2", "quat"]),
            _num(f"{prefix}_pose_lr", ["--pose-lr"], "1e-4"),
        ],
    )


ABINIT_GROUPS: list[Group] = [
    _g(
        "Checkpoint & seed",
        [
            _text("ab_load", ["--load"]),
            _text("ab_load_poses", ["--load-poses"]),
            _num("ab_seed", ["--seed"]),
            _flag("ab_verbose", ["-v"]),
        ],
    ),
    _g(
        "Dataset loading",
        [
            _text("ab_ind", ["--ind"]),
            _flag("ab_relion31", ["--relion31"]),
            _flag("ab_uninvert", ["--uninvert-data"]),
            _flag("ab_lazy", ["--lazy"]),
            _num("ab_max_threads", ["--max-threads"], "16"),
        ],
    ),
    _g(
        "Logging",
        [
            _num("ab_log_interval", ["--log-interval"], "10000"),
            _num("ab_checkpoint", ["--checkpoint"], "5"),
            _flag("ab_verbose_time", ["--verbose-time"]),
        ],
    ),
    _g(
        "Training parameters",
        [
            _num("ab_n", ["-n", "--num-epochs"], "30"),
            _num("ab_epochs_pose_search", ["--epochs-pose-search"]),
            _num("ab_n_imgs_pose_search", ["--n-imgs-pose-search"]),
            _num("ab_epochs_sgd", ["--epochs-sgd"]),
            _num("ab_pose_only_phase", ["--pose-only-phase"], "0"),
            _flag_false("ab_shuffle", ["--no-shuffle"]),
            _num("ab_num_workers", ["--num-workers"], "2"),
            _num("ab_shuffler_size", ["--shuffler-size"], "32768"),
            _flag("ab_multigpu", ["--multigpu"]),
            _no_amp("ab_use_amp"),
            _num("ab_batch_size_hps", ["--batch-size-hps"], "16"),
            _num("ab_batch_size_known_poses", ["--batch-size-known-poses"], "64"),
            _num("ab_batch_size_sgd", ["--batch-size-sgd"], "128"),
        ],
    ),
    _g(
        "Optimizers",
        [
            _num("ab_lr", ["--lr"], "1e-4"),
            _num("ab_lr_pose_table", ["--lr-pose-table"], "1e-3"),
            _num("ab_lr_conf_table", ["--lr-conf-table"], "1e-2"),
            _num("ab_lr_conf_encoder", ["--lr-conf-encoder"], "1e-4"),
            _num("ab_wd", ["--wd"], "0"),
            _select(
                "ab_hypervolume_optimizer_type",
                ["--hypervolume-optimizer-type"],
                ["adam"],
            ),
            _select(
                "ab_pose_table_optimizer_type",
                ["--pose-table-optimizer-type"],
                ["adam", "lbfgs"],
            ),
            _select(
                "ab_conf_table_optimizer_type",
                ["--conf-table-optimizer-type"],
                ["adam", "lbfgs"],
            ),
            _select(
                "ab_conf_encoder_optimizer_type",
                ["--conf-encoder-optimizer-type"],
                ["adam"],
            ),
        ],
    ),
    _g(
        "Masking",
        [
            _select(
                "ab_output_mask", ["--output-mask"], ["circ", "frequency_marching"]
            ),
            _num("ab_add_one_frequency_every", ["--add-one-frequency-every"], "100000"),
            _num("ab_n_frequencies_per_epoch", ["--n-frequencies-per-epoch"], "10"),
            _num("ab_max_freq", ["--max-freq"]),
            _num("ab_window_radius_gt_real", ["--window-radius-gt-real"], "0.85"),
        ],
    ),
    _g(
        "Losses",
        [
            _num("ab_beta_conf", ["--beta-conf"], "0"),
            _num("ab_trans_l1_regularizer", ["--trans-l1-regularizer"], "0"),
            _num("ab_l2_smoothness_regularizer", ["--l2-smoothness-regularizer"], "0"),
        ],
    ),
    _g(
        "Z / heterogeneity",
        [
            _flag("ab_variational_het", ["--variational-het"]),
            _num("ab_std_z_init", ["--std-z-init"], "0.1"),
            _flag("ab_use_conf_encoder", ["--use-conf-encoder"]),
            _num("ab_depth_cnn", ["--depth-cnn"], "5"),
            _num("ab_channels_cnn", ["--channels-cnn"], "32"),
            _num("ab_kernel_size_cnn", ["--kernel-size-cnn"], "3"),
            _num("ab_resolution_encoder", ["--resolution-encoder"]),
        ],
    ),
    _g(
        "Hypervolume",
        [
            _flag("ab_explicit_volume", ["--explicit-volume"]),
            _num("ab_layers", ["--layers"], "3"),
            _num("ab_dim", ["--dim"], "256"),
            _select("ab_pe_type", ["--pe-type"], ["gaussian"]),
            _num("ab_pe_dim", ["--pe-dim"], "64"),
            _num("ab_feat_sigma", ["--feat-sigma"], "0.5"),
            _select("ab_hypervolume_domain", ["--hypervolume-domain"], ["hartley"]),
            _select(
                "ab_pe_type_conf",
                ["--pe-type-conf"],
                ["", "geom"],
                labels=["(default None)", "geom"],
            ),
            _text("ab_initial_conf", ["--initial-conf"]),
        ],
    ),
    _g("Pretrain", [_num("ab_n_imgs_pretrain", ["--n-imgs-pretrain"], "10000")]),
    _g(
        "Pose search",
        [
            _num("ab_l_start", ["--l-start"], "12"),
            _num("ab_l_end", ["--l-end"], "32"),
            _num("ab_niter", ["--niter"], "4"),
            _num("ab_t_extent", ["--t-extent"], "20.0"),
            _num("ab_t_ngrid", ["--t-ngrid"], "7"),
            _num("ab_t_xshift", ["--t-xshift"], "0.0"),
            _num("ab_t_yshift", ["--t-yshift"], "0.0"),
            _flag(
                "ab_no_trans_search_at_pose_search",
                ["--no-trans-search-at-pose-search"],
            ),
            _num("ab_nkeptposes", ["--nkeptposes"], "8"),
            _num("ab_base_healpy", ["--base-healpy"], "2"),
            _flag("ab_no_trans", ["--no-trans"]),
        ],
    ),
    _g(
        "Normalization & analysis",
        [
            {"id": "ab_norm", "cli": ["--norm"], "w": "norm2"},
            _flag_false("ab_no_analysis", ["--no-analysis"]),
        ],
    ),
]

TRAIN_VAE_GROUPS: list[Group] = [
    _g(
        "I/O & logging",
        [
            _text("vae_load", ["--load"]),
            _flag_false("vae_no_analysis", ["--no-analysis"]),
            _num("vae_checkpoint", ["--checkpoint"], "1"),
            _num("vae_log_interval", ["--log-interval"], "1000"),
            _flag("vae_verbose", ["-v"]),
            _num("vae_seed", ["--seed"]),
            _num("vae_shuffle_seed", ["--shuffle-seed"]),
        ],
    ),
    _g(
        "Dataset loading",
        [
            _text("vae_ind", ["--ind"]),
            _flag("vae_uninvert", ["--uninvert-data"]),
            _flag_false("vae_window", ["--no-window"]),
            _num("vae_window_r", ["--window-r"], "0.85"),
            _flag("vae_lazy", ["--lazy"]),
            _num("vae_shuffler_size", ["--shuffler-size"], "0"),
            _num("vae_num_workers", ["--num-workers"], "0"),
            _num("vae_max_threads", ["--max-threads"], "16"),
        ],
    ),
    _g(
        "Tilt series parameters",
        [
            _num("vae_ntilts", ["--ntilts"], "10"),
            _flag("vae_random_tilts", ["--random-tilts"]),
            _num("vae_t_emb_dim", ["--t-emb-dim"], "64"),
            _num("vae_tlayers", ["--tlayers"], "3"),
            _num("vae_tdim", ["--tdim"], "1024"),
            _text("vae_dose_per_tilt", ["-d", "--dose-per-tilt"]),
            _num("vae_angle_per_tilt", ["-a", "--angle-per-tilt"], "3"),
        ],
    ),
    _g(
        "Training parameters",
        [
            _num("vae_n", ["-n", "--num-epochs"], "20"),
            _num("vae_b", ["-b", "--batch-size"], "16"),
            _num("vae_wd", ["--wd"], "0"),
            _num("vae_lr", ["--lr"], "1e-4"),
            _text("vae_beta", ["--beta"]),
            _text("vae_beta_control", ["--beta-control"]),
            {"id": "vae_norm", "cli": ["--norm"], "w": "norm2"},
            _no_amp("vae_use_amp"),
            _flag("vae_multigpu", ["--multigpu"]),
        ],
    ),
    _g(
        "Pose SGD",
        [
            _flag("vae_pose_sgd", ["--do-pose-sgd"]),
            _num("vae_pretrain", ["--pretrain"], "1"),
            _select("vae_emb_type", ["--emb-type"], ["s2s2", "quat"]),
            _num("vae_pose_lr", ["--pose-lr"], "3e-4"),
        ],
    ),
    _g(
        "Encoder Network",
        [
            _num("vae_enc_layers", ["--enc-layers"], "3"),
            _num("vae_enc_dim", ["--enc-dim"], "1024"),
            _select(
                "vae_encode_mode", ["--encode-mode"], ["conv", "resid", "mlp", "tilt"]
            ),
            _num("vae_enc_mask", ["--enc-mask"]),
            _flag("vae_use_real", ["--use-real"]),
        ],
    ),
    _g(
        "Decoder Network",
        [
            _num("vae_dec_layers", ["--dec-layers"], "3"),
            _num("vae_dec_dim", ["--dec-dim"], "1024"),
            _select(
                "vae_pe_type",
                ["--pe-type"],
                [
                    "geom_ft",
                    "geom_full",
                    "geom_lowf",
                    "geom_nohighf",
                    "linear_lowf",
                    "gaussian",
                    "none",
                ],
            ),
            _num("vae_feat_sigma", ["--feat-sigma"], "0.5"),
            _num("vae_pe_dim", ["--pe-dim"]),
            _select("vae_domain", ["--domain"], ["hartley", "fourier"]),
            _select("vae_activation", ["--activation"], ["relu", "leaky_relu"]),
        ],
    ),
]

TRAIN_NN_GROUPS: list[Group] = [
    _g(
        "I/O & logging",
        [
            _text("nn_load", ["--load"]),
            _num("nn_checkpoint", ["--checkpoint"], "1"),
            _num("nn_log_interval", ["--log-interval"], "1000"),
            _flag("nn_verbose", ["-v"]),
            _num("nn_seed", ["--seed"]),
            _num("nn_shuffle_seed", ["--shuffle-seed"]),
        ],
    ),
    _g(
        "Dataset loading",
        [
            _flag("nn_uninvert", ["--uninvert-data"]),
            _flag_false("nn_window", ["--no-window"]),
            _num("nn_window_r", ["--window-r"], "0.85"),
            _text("nn_ind", ["--ind"]),
            _flag("nn_lazy", ["--lazy"]),
            _num("nn_shuffler_size", ["--shuffler-size"], "0"),
        ],
    ),
    _g(
        "Training parameters",
        [
            _num("nn_n", ["-n", "--num-epochs"], "20"),
            _num("nn_b", ["-b", "--batch-size"], "16"),
            _num("nn_wd", ["--wd"], "0"),
            _num("nn_lr", ["--lr"], "1e-4"),
            {"id": "nn_norm", "cli": ["--norm"], "w": "norm2"},
            _no_amp("nn_use_amp"),
            _flag("nn_multigpu", ["--multigpu"]),
        ],
    ),
    _g(
        "Pose SGD",
        [
            _flag("nn_pose_sgd", ["--do-pose-sgd"]),
            _num("nn_pretrain", ["--pretrain"], "5"),
            _select("nn_emb_type", ["--emb-type"], ["s2s2", "quat"]),
            _num("nn_pose_lr", ["--pose-lr"], "1e-4"),
        ],
    ),
    _g(
        "Network Architecture",
        [
            _num("nn_layers", ["--layers"], "3"),
            _num("nn_dim", ["--dim"], "1024"),
            _num("nn_l_extent", ["--l-extent"], "0.5"),
            _select(
                "nn_pe_type",
                ["--pe-type"],
                [
                    "geom_ft",
                    "geom_full",
                    "geom_lowf",
                    "geom_nohighf",
                    "linear_lowf",
                    "gaussian",
                    "none",
                ],
            ),
            _num("nn_pe_dim", ["--pe-dim"]),
            _select("nn_domain", ["--domain"], ["hartley", "fourier"]),
            _select("nn_activation", ["--activation"], ["relu", "leaky_relu"]),
            _num("nn_feat_sigma", ["--feat-sigma"], "0.5"),
        ],
    ),
]

TRAIN_DEC_GROUPS: list[Group] = [
    _g(
        "I/O & logging",
        [
            _text("dec_load", ["--load"]),
            _flag_false("dec_no_analysis", ["--no-analysis"]),
            _num("dec_checkpoint", ["--checkpoint"], "1"),
            _num("dec_log_interval", ["--log-interval"], "1000"),
            _flag("dec_verbose", ["-v"]),
            _num("dec_seed", ["--seed"]),
            _num("dec_shuffle_seed", ["--shuffle-seed"]),
        ],
    ),
    _g(
        "Latent Variables",
        [
            _text("dec_load_z", ["--load-z"]),
            _num("dec_z_lr", ["--z-lr"], "1e-4"),
            _num("dec_pretrain_z", ["--pretrain-z"], "0"),
        ],
    ),
    _g(
        "Dataset loading",
        [
            _flag("dec_uninvert", ["--uninvert-data"]),
            _flag_false("dec_window", ["--no-window"]),
            _num("dec_window_r", ["--window-r"], "0.85"),
            _text("dec_ind", ["--ind"]),
            _flag("dec_lazy", ["--lazy"]),
            _num("dec_shuffler_size", ["--shuffler-size"], "0"),
        ],
    ),
    _g(
        "Training parameters",
        [
            _num("dec_n", ["-n", "--num-epochs"], "20"),
            _num("dec_b", ["-b", "--batch-size"], "8"),
            _num("dec_wd", ["--wd"], "0"),
            _num("dec_lr", ["--lr"], "1e-4"),
            {"id": "dec_norm", "cli": ["--norm"], "w": "norm2"},
            _no_amp("dec_use_amp"),
            _flag("dec_multigpu", ["--multigpu"]),
        ],
    ),
    _g(
        "Pose SGD",
        [
            _flag("dec_pose_sgd", ["--do-pose-sgd"]),
            _num("dec_pretrain_pose", ["--pretrain-pose"], "5"),
            _select("dec_emb_type", ["--emb-type"], ["s2s2", "quat"]),
            _num("dec_pose_lr", ["--pose-lr"], "1e-4"),
        ],
    ),
    _g(
        "Network Architecture",
        [
            _num("dec_layers", ["--layers"], "3"),
            _num("dec_dim", ["--dim"], "1024"),
            _num("dec_l_extent", ["--l-extent"], "0.5"),
            _select(
                "dec_pe_type",
                ["--pe-type"],
                [
                    "geom_ft",
                    "geom_full",
                    "geom_lowf",
                    "geom_nohighf",
                    "linear_lowf",
                    "gaussian",
                    "none",
                ],
            ),
            _num("dec_pe_dim", ["--pe-dim"]),
            _select("dec_domain", ["--domain"], ["hartley", "fourier"]),
            _select("dec_activation", ["--activation"], ["relu", "leaky_relu"]),
            _num("dec_feat_sigma", ["--feat-sigma"], "0.5"),
        ],
    ),
]

BACKPROJECT_VOXEL_GROUPS: list[Group] = [
    _g(
        "Dataset loading",
        [
            _flag("bpv_uninvert", ["--uninvert-data"]),
            _flag("bpv_lazy", ["--lazy"]),
            _text("bpv_ind", ["--ind"]),
            _num("bpv_first", ["--first"]),
        ],
    ),
    _g(
        "Backprojection parameters",
        [
            _flag_false("bpv_half_maps", ["--no-half-maps"]),
            _flag_false("bpv_fsc_vals", ["--no-fsc-vals"]),
            _num("bpv_batch_size", ["-b", "--batch-size"], "1000"),
            _select(
                "bpv_ctf_alg",
                ["--ctf-alg"],
                ["flip", "mul"],
                default="mul",
                help="CTF algorithm: phase flip (flip) or multiply (mul).",
            ),
            _num("bpv_reg_weight", ["--reg-weight"], "1.0"),
            _flag("bpv_output_sumcount", ["--output-sumcount"]),
            _num("bpv_log_interval", ["--log-interval"], "5000"),
        ],
    ),
    _g(
        "Tilt series parameters",
        [
            _flag("bpv_tilt", ["--tilt"]),
            _num("bpv_ntilts", ["--ntilts"], "10"),
            _flag("bpv_force_ntilts", ["--force-ntilts"]),
            _text("bpv_dose_per_tilt", ["-d", "--dose-per-tilt"]),
            _num("bpv_angle_per_tilt", ["-a", "--angle-per-tilt"], "3"),
        ],
    ),
]


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
            _num("ahet_l_start", ["--l-start"], "12"),
            _num("ahet_l_end", ["--l-end"], "32"),
            _num("ahet_niter", ["--niter"], "4"),
            _num("ahet_t_extent", ["--t-extent"], "10"),
            _num("ahet_t_ngrid", ["--t-ngrid"], "7"),
            _num("ahet_t_xshift", ["--t-xshift"], "0"),
            _num("ahet_t_yshift", ["--t-yshift"], "0"),
            _num("ahet_pretrain_ps", ["--pretrain"], "10000"),
            _num("ahet_ps_freq", ["--ps-freq"], "5"),
            _num("ahet_nkeptposes", ["--nkeptposes"], "8"),
            _num("ahet_base_healpy", ["--base-healpy"], "2"),
            _num("ahet_pose_model_update_freq", ["--pose-model-update-freq"]),
        ],
    )
    for i, g in enumerate(groups):
        if g["title"] == "Tilt series parameters":
            g["title"] = "Tilt series"
            g["args"].append(_flag("ahet_enc_only", ["--enc-only"]))
        elif g["title"] == "I/O & logging":
            g["args"].insert(1, _text("ahet_load_poses", ["--load-poses"]))
        elif g["title"] == "Training parameters":
            g["args"].extend(
                [
                    _text("ahet_equivariance", ["--equivariance"]),
                    _num("ahet_eq_start_it", ["--eq-start-it"]),
                    _num("ahet_eq_end_it", ["--eq-end-it"]),
                    _num("ahet_l_ramp_epochs", ["--l-ramp-epochs"], "0"),
                    _num("ahet_l_ramp_model", ["--l-ramp-model"], "0"),
                    _num("ahet_reset_model_every", ["--reset-model-every"]),
                    _num("ahet_reset_optim_every", ["--reset-optim-every"]),
                    _num(
                        "ahet_reset_optim_after_pretrain",
                        ["--reset-optim-after-pretrain"],
                    ),
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
            _text("ahom_tilt", ["--tilt"]),
            _num("ahom_tilt_deg", ["--tilt-deg"], "45"),
        ],
    )
    pose_search = _g(
        "Pose search parameters",
        [
            _num("ahom_l_start", ["--l-start"], "12"),
            _num("ahom_l_end", ["--l-end"], "32"),
            _num("ahom_niter", ["--niter"], "4"),
            _num("ahom_l_ramp_epochs", ["--l-ramp-epochs"], "25"),
            _flag("ahom_probabilistic", ["--probabilistic"]),
            _num("ahom_nkeptposes", ["--nkeptposes"], "8"),
            _num("ahom_base_healpy", ["--base-healpy"], "2"),
            _num("ahom_pose_model_update_freq", ["--pose-model-update-freq"]),
        ],
    )
    homo_train_extra: list[Arg] = [
        _num("ahom_t_extent", ["--t-extent"], "10"),
        _num("ahom_t_ngrid", ["--t-ngrid"], "7"),
        _num("ahom_t_xshift", ["--t-xshift"], "0"),
        _num("ahom_t_yshift", ["--t-yshift"], "0"),
        _flag("ahom_no_trans", ["--no-trans"]),
        _num("ahom_pretrain_train", ["--pretrain"], "10000"),
        _num("ahom_ps_freq", ["--ps-freq"], "5"),
    ]
    for i, g in enumerate(groups):
        if g["title"] == "I/O & logging":
            g["args"].insert(1, _text("ahom_load_poses", ["--load-poses"]))
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
            _num("ana_device", ["--device"]),
            _flag("ana_skip_vol", ["--skip-vol"]),
            _flag("ana_skip_umap", ["--skip-umap"]),
        ],
    ),
    _g(
        "Volumes to generate",
        [
            _num("ana_pc", ["--pc"], "2"),
            _num("ana_n_per_pc", ["--n-per-pc"], "10"),
            _num("ana_ksample", ["--ksample"], "20"),
        ],
    ),
    _g(
        "Volume post-processing",
        [
            _text("ana_apix", ["--Apix"]),
            _flag("ana_flip", ["--flip"]),
            _flag("ana_invert", ["--invert"]),
            _num("ana_downsample", ["-d", "--downsample"]),
            _text("ana_low_pass", ["--low-pass"]),
            _num("ana_crop", ["--crop"]),
            _num("ana_vol_start_index", ["--vol-start-index"], "1"),
        ],
    ),
]

ANALYZE_LANDSCAPE_GROUPS: list[Group] = [
    _g(
        "Run options",
        [
            _num("alsc_device", ["--device"]),
            _flag("alsc_multigpu", ["--multigpu"]),
            _flag("alsc_skip_umap", ["--skip-umap"]),
            _text("alsc_vol_ind", ["--vol-ind"]),
        ],
    ),
    _g(
        "Extra arguments for volume generation",
        [
            _num("alsc_sketch_size", ["-N", "--sketch-size"], "1000"),
            _num("alsc_apix", ["--Apix"], "1"),
            _flag("alsc_flip", ["--flip"]),
            _num("alsc_downsample", ["-d", "--downsample"], "128"),
            _flag("alsc_skip_vol", ["--skip-vol"]),
            _num("alsc_vol_start_index", ["--vol-start-index"], "1"),
        ],
    ),
    _g(
        "Extra arguments for mask generation",
        [
            _text("alsc_thresh", ["--thresh"]),
            _num("alsc_dilate", ["--dilate"], "5"),
            _num("alsc_cosine_edge", ["--cosine-edge"], "0"),
            _text("alsc_mask", ["--mask"]),
        ],
    ),
    _g(
        "Extra arguments for clustering",
        [
            _num("alsc_linkage", ["--linkage"], "average"),
            _num("alsc_n_clusters", ["-M"], "10"),
        ],
    ),
    _g(
        "Extra arguments for landscape visualization",
        [
            _num("alsc_pc_dim", ["--pc-dim"], "20"),
            _num("alsc_plot_dim", ["--plot-dim"], "5"),
        ],
    ),
]

ANALYZE_LANDSCAPE_FULL_GROUPS: list[Group] = [
    _g(
        "Run options",
        [
            _num("alfull_device", ["--device"]),
            _text("alfull_landscape_dir", ["--landscape-dir"]),
            _num("alfull_seed", ["--seed"], "0"),
        ],
    ),
    _g(
        "Volume generation arguments",
        [
            _num("alfull_training_volumes", ["-N", "--training-volumes"], "10000"),
            _flag("alfull_flip", ["--flip"]),
            _num("alfull_downsample", ["-d", "--downsample"], "128"),
            _flag("alfull_skip_vol", ["--skip-vol"]),
        ],
    ),
    _g(
        "Volume mapping arguments",
        [
            _num("alfull_batch_size", ["--batch-size"], "64"),
            _num("alfull_test_batch_size", ["--test-batch-size"], "1000"),
            _num("alfull_epochs", ["--epochs"], "200"),
            _num("alfull_lr", ["--lr"], "1e-4"),
            _num("alfull_dim", ["--dim"], "512"),
            _num("alfull_layers", ["--layers"], "3"),
        ],
    ),
    _g(
        "Volume PC clustering arguments",
        [
            _num("alfull_num_neighbors", ["--num-neighbors"], "50"),
            _num("alfull_resolution", ["--resolution"], "1.5"),
        ],
    ),
]

# Primary manuscript per reconstruction command (GitHub Pages nav link).
COMMAND_BUILDER_MANUSCRIPT_URLS: dict[str, str] = {
    "abinit": "https://www.nature.com/articles/s41592-025-02720-4",
    "train_dec": "https://www.nature.com/articles/s41592-025-02720-4",
    "train_vae": "https://www.nature.com/articles/s41592-020-01049-4",
    "train_nn": "https://www.nature.com/articles/s41592-020-01049-4",
    "abinit_het_old": (
        "https://openaccess.thecvf.com/content/ICCV2021/papers/"
        "Zhong_CryoDRGN2_Ab_Initio_Neural_Reconstruction_of_3D_Protein_Structures_From_"
        "ICCV_2021_paper.pdf"
    ),
    "abinit_homo_old": (
        "https://openaccess.thecvf.com/content/ICCV2021/papers/"
        "Zhong_CryoDRGN2_Ab_Initio_Neural_Reconstruction_of_3D_Protein_Structures_From_"
        "ICCV_2021_paper.pdf"
    ),
}

# Two-line nav label under the manuscript icon (GitHub Pages command builder).
COMMAND_BUILDER_MANUSCRIPT_LABELS: dict[str, str] = {
    "train_nn": "cryoDRGN1\nmanuscript",
    "train_vae": "cryoDRGN1\nmanuscript",
    "abinit_het_old": "cryoDRGN2\nmanuscript",
    "abinit_homo_old": "cryoDRGN2\nmanuscript",
    "abinit": "cryoDRGN-AI\nmanuscript",
    "train_dec": "cryoDRGN-AI\nmanuscript",
}

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

COMMAND_BUILDER_BATCH_SIZE_ARG_IDS: dict[str, list[str]] = {
    cmd: batch_size_arg_ids_for_command(cmd) for cmd in COMMAND_BUILDER_COMMAND_KEYS
}

_cli_help = load_cli_help_maps()
_cli_defaults = load_cli_defaults_maps()


def _attach_help(cmd: str, groups: list[Group]) -> None:
    attach_help_to_groups(
        _cli_help.get(cmd, {}),
        groups,
        defaults_map=_cli_defaults.get(cmd, {}),
    )


_attach_help("abinit", ABINIT_GROUPS)
_attach_help("abinit_het_old", ABINIT_HET_OLD_GROUPS)
_attach_help("abinit_homo_old", ABINIT_HOMO_OLD_GROUPS)
_attach_help("train_vae", TRAIN_VAE_GROUPS)
_attach_help("train_nn", TRAIN_NN_GROUPS)
_attach_help("train_dec", TRAIN_DEC_GROUPS)
_attach_help("backproject_voxel", BACKPROJECT_VOXEL_GROUPS)
_attach_help("analyze", ANALYZE_GROUPS)
_attach_help("analyze_landscape", ANALYZE_LANDSCAPE_GROUPS)
_attach_help("analyze_landscape_full", ANALYZE_LANDSCAPE_FULL_GROUPS)
attach_group_descriptions(COMMAND_BUILDER_SCHEMA)


def _required_field_titles(
    hm: dict[str, str],
    dm: dict[str, Any],
    pairs: dict[str, str | tuple[str, ...]],
) -> dict[str, str]:
    out: dict[str, str] = {}
    for elid, keys in pairs.items():
        seq: tuple[str, ...] = (keys,) if isinstance(keys, str) else keys
        for k in seq:
            t = resolved_help_for_flag(hm, dm, k)
            if t:
                out[elid] = t
                break
    return out


def _build_required_field_titles() -> dict[str, str]:
    r: dict[str, str] = {}
    r.update(
        _required_field_titles(
            _cli_help.get("abinit", {}),
            _cli_defaults.get("abinit", {}),
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
            _cli_defaults.get("abinit_het_old", {}),
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
            _cli_defaults.get("abinit_homo_old", {}),
            {
                "ahom_particles": "particles",
                "ahom_out": ("-o", "--outdir"),
            },
        ),
    )
    r.update(
        _required_field_titles(
            _cli_help.get("train_vae", {}),
            _cli_defaults.get("train_vae", {}),
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
            _cli_defaults.get("train_nn", {}),
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
            _cli_defaults.get("train_dec", {}),
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
            _cli_defaults.get("backproject_voxel", {}),
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
            _cli_defaults.get("analyze", {}),
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
            _cli_defaults.get("analyze_landscape", {}),
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
            _cli_defaults.get("analyze_landscape_full", {}),
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
        dm = _cli_defaults.get(cmd, {})
        if "--ctf" in hm:
            t = resolved_help_for_flag(hm, dm, "--ctf")
            if t:
                r[f"{prefix}_ctf"] = t
        if "--datadir" in hm:
            t = resolved_help_for_flag(hm, dm, "--datadir")
            if t:
                r[f"{prefix}_datadir"] = t
        if "--ind" in hm:
            t = resolved_help_for_flag(hm, dm, "--ind")
            if t:
                r[f"{prefix}_ind"] = t
    return r


COMMAND_BUILDER_REQUIRED_FIELD_TITLES: dict[str, str] = _build_required_field_titles()
