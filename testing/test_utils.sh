set -e
set -x

cryodrgn parse_ctf_star data/toy_projections.star -D 30 --Apix 1 -o test_ctf.pkl
cryodrgn_utils write_star data/toy_projections.mrcs --ctf test_ctf.pkl -o output/test.star
cryodrgn_utils write_star data/toy_projections.mrcs --ctf test_ctf.pkl -o output/test100.star --ind data/ind100.pkl
cryodrgn_utils write_star data/toy_projections.mrcs --ctf test_ctf.pkl -o output/test.star --poses data/toy_rot_trans.pkl

cryodrgn parse_pose_star output/test.star -D 30 -o output/test_pose.pkl
python diff_cryodrgn_pkl.py output/test_pose.pkl data/toy_rot_trans.pkl

# Test parsing from a .txt file
cryodrgn downsample data/toy_projections.mrcs -D 28 --chunk 80 -o data/toy_projections.mrcs
cryodrgn_utils write_star data/toy_projections.txt --ctf test_ctf.pkl -o output/test2.star
cryodrgn_utils write_star data/toy_projections.txt --ctf test_ctf.pkl --ind data/ind100.pkl -o output/test2_100.star

# Test copying micrograph coordinates
cryodrgn_utils write_star data/relion31.mrcs --ctf data/ctf1.pkl -o output/test3.star
#diff output/test3.star data/test3.star
cryodrgn_utils write_star data/relion31.mrcs --ctf data/ctf1.pkl -o output/test3_filtered.star --ind data/ind3.pkl
#diff output/test3_filtered.star data/test3_filtered.star
