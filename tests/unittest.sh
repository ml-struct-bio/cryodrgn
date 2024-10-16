# Script for running tests of cryoDRGN training and analysis methods outside of pytest
# NOTE: must be run within the folder containing this script:
#   $ cd cryodrgn/tests; sh unittest.sh

set -e
set -x

# first remove any output previously created using this script
rm -rf output/

# Test various file format loading options -- loss should be around 0.02
cryodrgn train_nn data/toy_projections.mrcs --poses data/toy_angles.pkl -o output/toy_recon -n 10 --no-amp
cryodrgn train_nn data/toy_projections.star --poses data/toy_angles.pkl -o output/toy_recon -n 10 --no-amp
cryodrgn train_nn data/toy_projections.txt --poses data/toy_angles.pkl -o output/toy_recon -n 10 --no-amp

# Test translations
cryodrgn_utils translate_mrcs data/toy_projections.mrcs data/toy_trans.pkl \
                              -o output/toy_projections.trans.mrcs --tscale -1
cryodrgn train_nn output/toy_projections.trans.mrcs --poses data/toy_rot_trans.pkl -o output/toy_recon -n 10 --no-amp
cryodrgn train_nn data/toy_projections.mrcs --poses data/toy_rot_zerotrans.pkl -o output/toy_recon -n 10 --no-amp

# Do pose SGD
cryodrgn train_nn data/toy_projections.mrcs --poses data/toy_rot_zerotrans.pkl \
                  -o output/toy_recon --do-pose-sgd --domain hartley

# Other decoder architectures
cryodrgn train_nn data/toy_projections.mrcs --poses data/toy_angles.pkl \
                  -o output/toy_recon --domain hartley -n 2 --no-amp
#cryodrgn train_nn data/toy_projections.mrcs --poses data/toy_angles.pkl \
#                  -o output/toy_recon --pe-type none -n 2 --no-amp
cryodrgn train_nn data/toy_projections.mrcs --poses data/toy_angles.pkl \
                  -o output/toy_recon --pe-type none --domain hartley -n 2 --no-amp

# Voxel-based backprojection
cryodrgn backproject_voxel data/hand.mrcs --poses data/hand_rot.pkl -o output/backproject.mrc
cryodrgn backproject_voxel data/sta_testing_bin8.star --ctf data/sta_ctf.pkl --poses data/sta_pose.pkl \
                           -o output/backproject_tilt.mrc --tilt --dose-per-tilt 3.0

# VAE
cryodrgn train_vae data/toy_projections.mrcs \
                   -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/toy_angles.pkl --zdim 10 -n 10
cryodrgn train_vae data/hand.mrcs -o output/hand_recon_vae --lr .0001 --seed 0 --poses data/hand_rot.pkl --zdim 10 -n 10
#cryodrgn train_vae data/hand.mrcs -o output/toy_recon_vae --lr .0001 --seed 0 --poses data/hand_rot.pkl \
#                                  --encode-mode conv --zdim 10 -n 10

# Test evaluation script
cryodrgn analyze output/toy_recon_vae 9
cryodrgn eval_vol output/toy_recon_vae/weights.pkl -c output/toy_recon_vae/config.yaml \
                  -z 0 0 0 0 0 0 0 0 0 0 -o output/toy_recon_vae/vol.mrc
cryodrgn eval_images data/toy_projections.mrcs output/toy_recon_vae/weights.pkl -c output/toy_recon_vae/config.yaml \
                     -o output/toy_recon_vae/losses.pkl --out-z output/toy_recon_vae/z_eval.pkl \
                     --poses data/toy_angles.pkl
cryodrgn pc_traversal output/toy_recon_vae/z.pkl -o output/toy_recon_vae/pc_traversal
cryodrgn graph_traversal output/toy_recon_vae/z.pkl --anchors 0 10 100 \
                         -o output/toy_recon_vae/graph_traversal/path.txt \
                         --outind output/toy_recon_vae/graph_traversal/z.path.txt

# Test apex.amp
cryodrgn train_nn data/hand.mrcs --poses data/hand_rot.pkl -o output/hand_recon -b 8 --no-amp
cryodrgn train_nn data/hand.mrcs --poses data/hand_rot.pkl -o output/hand_recon -b 16

# CTF testing
cryodrgn parse_ctf_csparc data/cryosparc_P12_J24_001_particles.cs -o test_ctf.pkl
cryodrgn parse_ctf_star data/toy_projections.star -D 30 --Apix 1 -o test_ctf.pkl
cryodrgn train_vae data/toy_projections.mrcs -o output/toy_recon_vae \
                   --lr .0001 --seed 0 --poses data/toy_angles.pkl --ctf test_ctf.pkl --zdim 10

set +x
echo ">>>>>   All unittest.sh tests passed!   <<<<<"
