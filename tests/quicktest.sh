#!/usr/bin/env sh
set -e
set -x
cryodrgn train_vae  data/hand.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/hand_rot.pkl --zdim 10 --pe-type gaussian
cryodrgn train_vae  data/hand.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/hand_rot.pkl --zdim 10 --pe-type gaussian --load latest -n 30
