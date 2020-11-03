set -e 
set -x

cryodrgn parse_ctf_star data/toy_projections.star -D 30 --Apix 1 -o test_ctf.pkl
python ../utils/write_starfile.py data/toy_projections.mrcs test_ctf.pkl -o output/test.star
python ../utils/write_starfile.py data/toy_projections.mrcs test_ctf.pkl -o output/test100.star --ind data/ind100.pkl 
python ../utils/write_starfile.py data/toy_projections.mrcs test_ctf.pkl -o output/test.star --poses data/toy_rot_trans.pkl

cryodrgn parse_pose_star output/test.star -D 30 -o output/test_pose.pkl
#diff output/test_pose.pkl data/toy_rot_trans.pkl
