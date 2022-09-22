set -e 
set -x

cryodrgn parse_ctf_star data/toy_projections.star -D 30 --Apix 1 -o test_ctf.pkl
cryodrgn write_starfile data/toy_projections.mrcs test_ctf.pkl -o output/test.star
cryodrgn write_starfile data/toy_projections.mrcs test_ctf.pkl -o output/test100.star --ind data/ind100.pkl 
cryodrgn write_starfile data/toy_projections.mrcs test_ctf.pkl -o output/test.star --poses data/toy_rot_trans.pkl

cryodrgn parse_pose_star output/test.star -D 30 -o output/test_pose.pkl
python diff_cryodrgn_pkl.py output/test_pose.pkl data/toy_rot_trans.pkl

# Test parsing from a .txt file 
cryodrgn downsample data/toy_projections.mrcs -D 28 --chunk 80 -o data/toy_projections.mrcs
cryodrgn write_starfile data/toy_projections.txt test_ctf.pkl -o output/test2.star
cryodrgn write_starfile data/toy_projections.txt test_ctf.pkl --ind data/ind100.pkl -o output/test2_100.star

# Test copying micrograph coordinates
cryodrgn write_starfile data/relion31.mrcs data/ctf1.pkl -o output/test3.star --ref-star data/FinalRefinement-OriginalParticles-PfCRT.star --keep-micrograph 
#diff output/test3.star data/test3.star
cryodrgn write_starfile data/relion31.mrcs data/ctf1.pkl -o output/test3_filtered.star --ref-star data/FinalRefinement-OriginalParticles-PfCRT.star --keep-micrograph --ind data/ind3.pkl 
#diff output/test3_filtered.star data/test3_filtered.star
