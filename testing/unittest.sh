set -e
set -x

python ../vae_rot.py  data/toy_projections.mrcs -o output/toy_recon_vae --lr .0001 --seed 0
python ../vae_rot.py  data/toy_projections.mrcs -o output/toy_recon_vae --lr .0001 --no-trans --seed 0

python ../vae_tilt.py  data/toy_projections.mrcs data/toy_projections.mrcs --tilt 45 -o output/toy_recon_vae --lr .0001 --seed 0
python ../vae_tilt.py  data/toy_projections.mrcs data/toy_projections.mrcs --tilt 45 -o output/toy_recon_vae --lr .0001 --no-trans --seed 0

python ../bnb_rot.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0 --l-start 10 --l-end 14 -b 10
python ../bnb_rot.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0 --no-trans --l-start 10 --l-end 14 -b 10
python ../bnb_rot.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0 --tilt data/toy_projections.mrcs --tilt-deg 45 --l-start 10 --l-end 14 -b 10
python ../bnb_rot.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0 --no-trans --tilt data/toy_projections.mrcs --tilt-deg 45 --l-start 10 --l-end 14 -b 10

python ../bnb_het.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0 -b 10 -n 1 --l-start 10 --l-end 14
python ../bnb_het.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0 -b 10 --tilt data/toy_projections.mrcs --tilt-deg 45 --encode-mode tilt -n 1 --l-start 10 --l-end 14
python ../bnb_het.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0 -b 10 --tilt data/toy_projections.mrcs --tilt-deg 45 --encode-mode tilt --rotate --enc-only -n 1 --l-start 10 --l-end 14
python ../bnb_het.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0 -b 10 --tilt data/toy_projections.mrcs --tilt-deg 45 --encode-mode tilt --enc-only -n 1 --l-start 10 --l-end 14

python ../vae_priors.py  data/toy_projections.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --priors data/toy_angles.pkl --no-trans
python ../vae_priors.py  data/toy_projections.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --priors data/toy_angles.pkl --no-trans --pretrain 2

python ../backproject_nn_hartley.py data/toy_projections.mrcs data/toy_angles.pkl --l-extent 3 -o output/toy_recon
python ../backproject_nn_hartley.py data/toy_projections.mrcs data/toy_angles.pkl --l-extent 3 -o output/toy_recon --trans data/tilt_series/trans.zero.pkl 

python ../backproject_k_enc2.py data/toy_projections.mrcs --poses data/toy_angles.pkl -o output/toy_recon
python ../backproject_k_enc2.py data/toy_projections.mrcs --poses data/toy_angles.pkl data/toy_trans.zero.pkl -o output/toy_recon
python ../backproject_k_enc2.py data/toy_projections.mrcs --poses data/toy_angles.pkl data/toy_trans.zero.pkl -o output/toy_recon --do-pose-sgd
python ../backproject_k_enc2.py data/toy_projections.mrcs --poses data/toy_angles.pkl -o output/toy_recon --domain hartley
python ../backproject_k_enc2.py data/toy_projections.mrcs --poses data/toy_angles.pkl -o output/toy_recon --pe-type none
python ../backproject_k_enc2.py data/toy_projections.mrcs --poses data/toy_angles.pkl -o output/toy_recon --pe-type none --domain hartley

python ../backproject.py data/hand.mrcs --poses data/hand_rot.pkl -o output/backproject.mrc
python ../backproject.py data/hand.mrcs --poses data/hand_rot.pkl -o output/backproject_tilt.mrc --tilt data/hand_tilt.mrcs

python ../vae_het.py  data/toy_projections.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/toy_angles.pkl
python ../vae_het.py  data/hand.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/hand_rot.pkl 
python ../vae_het.py  data/hand.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/hand_rot.pkl --encode-mode conv
