#!/usr/bin/env sh
set -e
set -x

# cryodrgn
cryodrgn train_vae data/hand.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 \
                                  --poses data/hand_rot_trans.pkl --ctf data/test_ctf.100.pkl --zdim 10
cryodrgn analyze output/toy_recon_vae

# drgnai-fixed
cryodrgn setup output/toy_recon_vae_2 --particles data/hand.mrcs --poses data/hand_rot_trans.pkl \
                                      --zdim 10 --model cryodrgn-ai --pose-estimation fixed \
                                      --cfg hidden_dim=1024 --ctf data/test_ctf.100.pkl
cryodrgn train output/toy_recon_vae_2


# cryodrgn abinit
cryodrgn abinit_het data/hand.mrcs -o output/toy_recon_vae_ab --lr .0001 --seed 0 \
                                   --ctf data/test_ctf.100.pkl --zdim 10
cryodrgn analyze output/toy_recon_vae_ab

# drgnai-abinit
cryodrgn setup output/toy_recon_vae_ab2 --particles data/hand.mrcs --zdim 10 --model cryodrgn-ai \
                                        --cfg hidden_dim=1024 t_ngrid=3 n_imgs_pose_search=100 \
                                        --ctf data/test_ctf.100.pkl
cryodrgn train output/toy_recon_vae_ab2
