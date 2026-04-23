"""Structured optional-arg groups for the dashboard command builder.

Mirrors ``add_argument_group`` titles and flags from:
``cryodrgn.commands.abinit``, ``train_vae``, ``train_nn``, ``train_dec``.
``--ctf`` and ``--datadir`` are rendered in the template (required row and
top-row dataset column), not listed here. Other required fields are template-only.
"""

from __future__ import annotations

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


def _g(title: str, args: list[Arg]) -> Group:
    return {"title": title, "args": args}


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

COMMAND_BUILDER_SCHEMA: Schema = {
    "abinit": ABINIT_GROUPS,
    "train_vae": TRAIN_VAE_GROUPS,
    "train_nn": TRAIN_NN_GROUPS,
    "train_dec": TRAIN_DEC_GROUPS,
}

_cli_help = load_cli_help_maps()
attach_help_to_groups(_cli_help.get("abinit", {}), ABINIT_GROUPS)
attach_help_to_groups(_cli_help.get("train_vae", {}), TRAIN_VAE_GROUPS)
attach_help_to_groups(_cli_help.get("train_nn", {}), TRAIN_NN_GROUPS)
attach_help_to_groups(_cli_help.get("train_dec", {}), TRAIN_DEC_GROUPS)


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
    for prefix, cmd in (
        ("ab", "abinit"),
        ("vae", "train_vae"),
        ("nn", "train_nn"),
        ("dec", "train_dec"),
    ):
        hm = _cli_help.get(cmd, {})
        if "--ctf" in hm:
            r[f"{prefix}_ctf"] = hm["--ctf"]
        if "--datadir" in hm:
            r[f"{prefix}_datadir"] = hm["--datadir"]
    return r


COMMAND_BUILDER_REQUIRED_FIELD_TITLES: dict[str, str] = _build_required_field_titles()
