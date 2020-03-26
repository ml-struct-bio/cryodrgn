set -e
set -x

# Test various file format loading options -- loss should be around 0.02
python ../backproject_nn.py data/toy_projections.mrcs --poses data/toy_angles.pkl -o output/toy_recon
python ../backproject_nn.py data/toy_projections.star --poses data/toy_angles.pkl -o output/toy_recon
python ../backproject_nn.py data/toy_projections.txt --poses data/toy_angles.pkl -o output/toy_recon

# Test translations
python ../utils/translate_stack.py data/toy_projections.mrcs data/toy_trans.pkl -o output/toy_projections.trans.mrcs --tscale -1
python ../backproject_nn.py output/toy_projections.trans.mrcs --poses data/toy_poses_wtrans.pkl -o output/toy_recon
python ../backproject_nn.py data/toy_projections.mrcs --poses data/toy_rot_zerotrans.pkl -o output/toy_recon

# Do pose SGD
python ../backproject_nn.py data/toy_projections.mrcs --poses data/toy_rot_zerotrans.pkl -o output/toy_recon --do-pose-sgd

# Other decoder architectures
python ../backproject_nn.py data/toy_projections.mrcs --poses data/toy_angles.pkl -o output/toy_recon --domain hartley
python ../backproject_nn.py data/toy_projections.mrcs --poses data/toy_angles.pkl -o output/toy_recon --pe-type none
python ../backproject_nn.py data/toy_projections.mrcs --poses data/toy_angles.pkl -o output/toy_recon --pe-type none --domain hartley

# Voxel-based backprojection
python ../backproject_voxel.py data/hand.mrcs --poses data/hand_rot.pkl -o output/backproject.mrc
python ../backproject_voxel.py data/hand.mrcs --poses data/hand_rot.pkl -o output/backproject_tilt.mrc --tilt data/hand_tilt.mrcs

# VAE
python ../vae_het.py  data/toy_projections.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/toy_angles.pkl --zdim 10
python ../vae_het.py  data/hand.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/hand_rot.pkl --zdim 10 
python ../vae_het.py  data/hand.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/hand_rot.pkl --encode-mode conv --zdim 10

# Test apex.amp
python ../backproject_nn.py data/hand.mrcs --poses data/hand_rot.pkl -o output/hand_recon -b 8
python ../backproject_nn.py data/hand.mrcs --poses data/hand_rot.pkl -o output/hand_recon --amp -b 8

# CTF testing
python ../utils/parse_ctf_csparc.py data/cryosparc_P12_J24_001_particles.cs -o test_ctf.pkl
python ../utils/parse_ctf_star.py data/toy_projections.star -N 1000 --Apix 1 -o test_ctf.pkl
python ../vae_het.py  data/toy_projections.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/toy_angles.pkl --ctf test_ctf.pkl --zdim 10
