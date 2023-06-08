import time

import numpy as np
import torch

from cryodrgn import dataset, lattice, models, pose_search, utils

use_cuda = torch.cuda.is_available()
device = torch.device("cuda" if use_cuda else "cpu")
print("Use cuda {}".format(use_cuda))
# if use_cuda:
#     torch.set_default_tensor_type(torch.cuda.FloatTensor)  # type: ignore

basedir = "datasets/ribo_real_128"
data = dataset.ImageDataset(
    f"{basedir}/particles.128.phaseflip.1000.mrcs",
    lazy=False,
    window=False,
    keepreal=True,
)


S = 0
D = data.D
assert D % 2 == 1

lat = lattice.Lattice(D)
pose = utils.load_pkl(f"{basedir}/gt.pose.pkl")
pose_rot, pose_trans = pose
pose_rot = torch.tensor(pose_rot)
pose_trans = torch.tensor(pose_trans.astype(np.float32) * (D - 1))


def load_model(path):
    ckpt = torch.load(path)
    model = models.get_decoder(3, D, 3, 256, "fourier", "geom_lowf")
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()
    if use_cuda:
        model.cuda()
    return model


model = load_model(f"{basedir}/trained_models/weights_gt_poses.pkl")
print(f"Device: {next(model.parameters()).device}")


def do_pose_search(
    images, model, base_healpy=2, nkeptposes=8, Lmin=12, Lmax=24, niter=5, **kwargs
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

    return ps.opt_theta_trans(images)


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


def eval_pose_search(data, model, B=512, label="", **kwargs):
    tic = time.perf_counter()
    res = []
    for chunk in torch.from_numpy(data.particles[S : S + B]).split(8):
        res.append(do_pose_search(chunk, model, **kwargs))
    delta = time.perf_counter() - tic
    batch_rot = pose_rot[S : S + B]
    batch_trans = pose_trans[S : S + B]

    rot_hat, trans_hat, _ = [torch.cat(x) for x in zip(*res)]

    print(
        f"{label:20s}| "
        f"Rot MedSE= {medse(rot_hat, batch_rot):.5f} "
        f"Rot MSE= {mse(rot_hat, batch_rot):.5f} "
        f"Trans MedSE= {medse(trans_hat, batch_trans):.5f} "
        f"Trans MSE= {mse(trans_hat, batch_trans):.5f} "
        f"Trans offset= {trans_offset(trans_hat, batch_trans)} "
        f"time= {delta:.2f} s"
    )


print("=" * 80)

tic = time.perf_counter()

eval_pose_search(
    data,
    model,
    label="base",
)

eval_pose_search(
    data,
    model,
    label="L= [12, 48]",
    Lmin=12,
    Lmax=48,
)

eval_pose_search(
    data,
    model,
    label="L= [12, 48], 7 iters",
    Lmin=12,
    Lmax=48,
    niter=7,
)

eval_pose_search(
    data,
    model,
    label="L= [12, 48], 7 iters, half",
    Lmin=12,
    Lmax=48,
    niter=7,
    half_precision=True,
)

for t_ngrid in (3, 14):
    eval_pose_search(
        data,
        model,
        label=f"t_ngrid= {t_ngrid}",
        t_ngrid=t_ngrid,
    )

for niter in (2, 5, 7):
    eval_pose_search(
        data,
        model,
        label=f"niter= {niter}",
        niter=niter,
    )

for nkp in (1, 2, 4, 24):
    eval_pose_search(
        data,
        model,
        label=f"nkp= {nkp:2d}",
        nkeptposes=nkp,
    )

eval_pose_search(
    data,
    model,
    label="healpy= 3",
    base_healpy=3,
)

# for bhp in (1, 2, 3):
#     for nkp in (1, 4, 12, 24):
#         eval_pose_search(data, model,
#                         label=f"noisy nkp= {nkp:2d} bhp= {bhp}",
#                         base_healpy=bhp,
#                         nkeptposes=nkp)
#     print('-' * 80)

print(f"Finished in {time.perf_counter() - tic} s ")

# import cProfile
# pr = cProfile.Profile()
# pr.enable()
# eval_pose_search(data, model, nkeptposes=24)
# pr.disable()
# pr.print_stats('cumtime')
