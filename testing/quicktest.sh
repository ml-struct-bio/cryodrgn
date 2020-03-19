set -e
set -x

python $CDRGN_SRC/vae_het.py  $CDRGN_SRC/testing/data/hand.mrcs -o /tmp/toy_recon_vae --lr .0001 --seed 0 --poses $CDRGN_SRC/testing/data/hand_rot.pkl --zdim 10 
