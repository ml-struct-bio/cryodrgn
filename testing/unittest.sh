set -e
set -x

python ../vae_rot.py  data/toy_projections.mrcs -o output/toy_recon_vae --lr .0001
python ../vae_rot.py  data/toy_projections.mrcs -o output/toy_recon_vae --lr .0001 --no-trans
python ../vae_tilt.py  data/toy_projections.mrcs data/toy_projections.mrcs --tilt 45 -o output/toy_recon_vae --lr .0001 
python ../vae_tilt.py  data/toy_projections.mrcs data/toy_projections.mrcs --tilt 45 -o output/toy_recon_vae --lr .0001 --no-trans
