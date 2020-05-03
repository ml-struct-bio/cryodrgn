import cryodrgn
from cryodrgn import mrc, dataset, models, lattice, fft, utils, pose_search
import matplotlib.pyplot as plt
import seaborn as sns
import math
import time
import torch
import numpy as np
from scipy.spatial.transform import Rotation as RR

use_cuda = torch.cuda.is_available()
device = torch.device('cuda' if use_cuda else 'cpu')
print('Use cuda {}'.format(use_cuda))
if use_cuda:
    torch.set_default_tensor_type(torch.cuda.FloatTensor)

basedir = "datasets/ribo_syn_64"
data = dataset.MRCData(f'{basedir}/projections.1k.mrcs', window=False, keepreal=True)
data_noisy = dataset.MRCData(f'{basedir}/noise_0.1/projections.1k.mrcs', window=False, keepreal=True)

D = data.D
lat = lattice.Lattice(D)

pose = utils.load_pkl(f'{basedir}/pose.pkl')
pose_rot, pose_trans = pose
pose_rot = torch.tensor(pose_rot)
pose_trans = torch.tensor(pose_trans.astype(np.float32) * 64)

def load_model(path):
    ckpt = torch.load(path)
    model = models.get_decoder(3, 65, 3, 256, 'fourier', 'geom_lowf')
    model.load_state_dict(ckpt['model_state_dict'])
    model.eval()
    if use_cuda:
        model.cuda()
    return model

model = load_model(f'{basedir}/trained_model/weights.pkl')
model_noisy = load_model(f'{basedir}/trained_model_noise/weights.pkl')
print(f"Device: {next(model.parameters()).device}")

def do_pose_search(images, model, nkeptposes=24, Lmin=12, Lmax=24, niter=5, **kwargs):
    device = next(model.parameters()).device
    images = torch.from_numpy(images).to(device)
    ps = pose_search.PoseSearch(
        model=model,
        lattice=lat,
        Lmin=Lmin,
        Lmax=Lmax,
        nkeptposes=nkeptposes,
        **kwargs,
    )

    return ps.opt_theta_trans(images, niter=niter)

def mse(x, y):
    B = x.shape[0]
    return (x - y).pow(2).view(B, -1).sum(-1).mean()

def medse(x, y):
    B = x.shape[0]
    return (x - y).pow(2).view(B, -1).sum(-1).median()

def eval_pose_search(data, model, B=64, label="", **kwargs):
    rot_hat, trans_hat = do_pose_search(data.particles[:B], model, **kwargs)
    print(f"{label} "
          f"Rot MedSE= {medse(rot_hat, pose_rot[:B]):.4f} "
          f"Rot MSE= {mse(rot_hat, pose_rot[:B]):.4f} "
          f"Trans MedSE= {medse(trans_hat, pose_trans[:B]):.4f} "
          f"Trans MSE= {mse(trans_hat, pose_trans[:B]):.4f}")

print("=============================================")
# for nkp in (1, 4, 12):
#     eval_pose_search(data, model, label=f"nkp= {nkp}", nkeptposes=nkp)

# for nkp in (1, 4, 12, 24):
#     eval_pose_search(data_noisy, model_noisy,
#                      label=f"noisy nkp= {nkp}",
#                      nkeptposes=nkp)


tic = time.perf_counter()

for bhp in (1, 2, 3):
    for nkp in (1, 4, 12, 24):
        eval_pose_search(data_noisy, model_noisy,
                        B=8,
                        label=f"noisy nkp= {nkp} bhp= {bhp}",
                        base_healpy=bhp,
                        nkeptposes=nkp)

print(f"Finished in {time.perf_counter() - tic} s ")

# import cProfile
# pr = cProfile.Profile()
# pr.enable()
# eval_pose_search(data, model, nkeptposes=24)
# pr.disable()
# pr.print_stats('cumtime')
