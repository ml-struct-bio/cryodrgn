import argparse
import os
import time

import numpy as np
import torch

from cryodrgn import dataset, lattice, models, pose_search, utils

use_cuda = torch.cuda.is_available()
device = torch.device("cuda" if use_cuda else "cpu")
print("Use cuda {}".format(use_cuda))
# if use_cuda:
#     torch.set_default_tensor_type(torch.cuda.FloatTensor)  # type: ignore


def load_model(path, D):
    print(f"Loading model from {path}")
    ckpt = torch.load(path)
    model = models.get_decoder(3, D, 3, 256, "hartley", "geom_lowf")
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()
    if use_cuda:
        model.cuda()
    return model


def get_poses(path, D):
    print(f"Loading poses from {path}")
    pose = utils.load_pkl(path)
    pose_rot, pose_trans = pose
    pose_rot = torch.tensor(pose_rot)
    pose_trans = torch.tensor(pose_trans.astype(np.float32) * (D - 1))
    return pose_rot, pose_trans


def mse(x, y):
    B = x.shape[0]
    errors = (x - y).pow(2).view(B, -1).sum(-1)
    # print('mse', errors)
    return errors.mean()


def medse(x, y):
    B = x.shape[0]
    return (x - y).pow(2).view(B, -1).sum(-1).median()


def trans_offset(x, y):
    return (x - y).view(-1, 2).mean(0).cpu().numpy()


def run(args):
    GPU_BATCH = 4

    print(f"Loading particle images from {args.particles}")
    data = dataset.ImageDataset(
        args.particles, window=False, keepreal=False, invert_data=True
    )

    LB = ""
    D = data.D
    assert D % 2 == 1

    lat = lattice.Lattice(D)

    model = load_model(args.model, D)
    pose_rot, pose_trans = get_poses(args.poses, D)

    def do_pose_search(
        images,
        base_healpy=2,
        nkeptposes=8,
        Lmin=12,
        Lmax=24,
        niter=5,
        **kwargs,
    ):
        device = next(model.parameters()).device
        images = images.to(device)
        ps = pose_search.PoseSearch(
            model=model,
            lattice=lat,
            Lmin=Lmin,
            Lmax=Lmax,
            nkeptposes=nkeptposes,
            t_extent=10,
            base_healpy=base_healpy,
            niter=niter,
            **kwargs,
        )

        ret = ps.opt_theta_trans(images)
        return ret

    def eval_pose_search(B=256, S=0, label="", **kwargs):
        tic = time.perf_counter()
        res = []
        particles, _, _ = data[S : S + B]
        for chunk in particles.split(GPU_BATCH):
            res.append(do_pose_search(chunk, **kwargs))
        delta = time.perf_counter() - tic
        batch_rot = pose_rot[S : S + B]
        batch_trans = pose_trans[S : S + B]

        rot_hat, trans_hat, _ = [torch.cat(x) for x in zip(*res)]

        # print(rot_hat)
        # print(batch_rot)

        # print(pose_trans[:10])
        # print(trans_hat[:10])

        print(
            f"{label:20s}| "
            f"Rot MedSE= {medse(rot_hat, batch_rot):.5f} "
            f"Rot MSE= {mse(rot_hat, batch_rot):.5f} "
            f"Trans MedSE= {medse(trans_hat, batch_trans):.5f} "
            f"Trans MSE= {mse(trans_hat, batch_trans):.5f} "
            f"Trans offset= {trans_offset(trans_hat, batch_trans)} "
            f"time= {delta:.2f} s"
        )

    print(f"Device: {next(model.parameters()).device}")

    print("=" * 80)

    tic = time.perf_counter()

    eval_pose_search(
        label=f"{LB}base",
    )

    # eval_pose_search(
    #     label=f"{LB}L= [12, 48]",
    #     Lmin=12,
    #     Lmax=48,
    # )

    # eval_pose_search(
    #     label=f"{LB}L= [12, 48], 7 iters",
    #     Lmin=12,
    #     Lmax=48,
    #     niter=7,
    # )

    # eval_pose_search(
    #     label=f"{LB}L= [12, 48], 7 iters, half",
    #     Lmin=12,
    #     Lmax=48,
    #     niter=7,
    #     half_precision=True,
    # )

    # for t_ngrid in (3, 14):
    #     eval_pose_search(
    #         label=f"{LB}t_ngrid= {t_ngrid}",
    #         t_ngrid=t_ngrid,
    #     )

    # for niter in (2, 5, 7):
    #     eval_pose_search(
    #         label=f"{LB}niter= {niter}",
    #         niter=niter,
    #     )

    # for nkp in (1, 2, 4, 24):
    #     eval_pose_search(
    #         label=f"{LB}nkp= {nkp:2d}",
    #         nkeptposes=nkp,
    #     )

    # eval_pose_search(
    #     label=f"{LB}healpy= 3",
    #     base_healpy=3,
    # )

    print(f"Finished in {time.perf_counter() - tic} s ")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--particles",
        default=f"{os.path.dirname(__file__)}/../datasets/rag12_128/particles.128.phaseflip.1k.mrcs",
    )
    parser.add_argument(
        "--model",
        default=f"{os.path.dirname(__file__)}/../datasets/rag12_128/good_model/weights.50.pkl",
    )
    parser.add_argument(
        "--poses",
        default=f"{os.path.dirname(__file__)}/../datasets/rag12_128/good_model/pose.v2.pkl",
    )

    parser.add_argument("--B", type=int, default=512)
    args = parser.parse_args()

    run(args)
