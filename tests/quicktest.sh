#!/usr/bin/env sh
set -e
set -x

# cryodrgn
#cryodrgn train_vae  data/hand.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/hand_rot.pkl --zdim 10 --pe-type gaussian
cryodrgn train_vae  data/hand.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/hand_rot_trans.pkl --ctf data/test_ctf.100.pkl --zdim 10
cryodrgn analyze output/toy_recon_vae

# drgnai-fixed
cryodrgn setup output/toy_recon_vae_2 --particles data/hand.mrcs --poses data/hand_rot_trans.pkl --z-dim 10 --model cryodrgn-ai --pose-estimation fixed --cfg n_imgs_pose_search=50 --ctf data/test_ctf.100.pkl 
cryodrgn train output/toy_recon_vae_2
