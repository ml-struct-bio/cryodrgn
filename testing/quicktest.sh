set -e
set -x

python ../vae_het.py  data/hand.mrcs -o /tmp/toy_recon_vae --lr .0001 --seed 0 --poses data/hand_rot.pkl --zdim 10 
