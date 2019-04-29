set -e
set -x

python ../vae_rot.py  data/toy_projections.mrcs -o output/toy_recon_vae --lr .0001 --seed 0
python ../vae_rot.py  data/toy_projections.mrcs -o output/toy_recon_vae --lr .0001 --no-trans --seed 0
python ../vae_tilt.py  data/toy_projections.mrcs data/toy_projections.mrcs --tilt 45 -o output/toy_recon_vae --lr .0001 --seed 0
python ../vae_tilt.py  data/toy_projections.mrcs data/toy_projections.mrcs --tilt 45 -o output/toy_recon_vae --lr .0001 --no-trans --seed 0
python ../bnb_rot.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0
python ../bnb_rot.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0 --no-trans
python ../bnb_rot.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0 --tilt data/toy_projections.mrcs --tilt-deg 45
python ../bnb_rot.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0 --no-trans --tilt data/toy_projections.mrcs --tilt-deg 45
python ../bnb_het.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0 -b 10 -n 1
python ../bnb_het.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0 -b 10 --tilt data/toy_projections.mrcs --tilt-deg 45 --encode-mode tilt -n 1
python ../bnb_het.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0 -b 10 --tilt data/toy_projections.mrcs --tilt-deg 45 --encode-mode tilt --rotate --enc-only -n 1
python ../bnb_het.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0 -b 10 --tilt data/toy_projections.mrcs --tilt-deg 45 --encode-mode tilt --enc-only -n 1
