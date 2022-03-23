set -e
set -x

# Test 1: Pase CTF from starfile with missing parameter
cryodrgn parse_ctf_star data/FinalRefinement-OriginalParticles-PfCRT.star -w .1 -o output/ctf1.pkl --png output/ctf1.png -D 300 --Apix 1.035
python diff_cryodrgn_pkl.py output/ctf1.pkl data/ctf1.pkl

# Test 2: Parse CTF from cs file
cryodrgn parse_ctf_csparc data/cryosparc_P12_J24_001_particles.cs -o output/ctf2.pkl --png output/ctf2.png
python diff_cryodrgn_pkl.py output/ctf2.pkl data/ctf2.pkl

# Test 3: Parse poses from starfile
cryodrgn parse_pose_star data/FinalRefinement-OriginalParticles-PfCRT.star -o output/pose.star.pkl -D 300
python diff_cryodrgn_pkl.py output/pose.star.pkl data/pose.star.pkl

# Test 4: Parse poses from .cs file
cryodrgn parse_pose_csparc data/cryosparc_P12_J24_001_particles.cs -D 180 -o output/pose.cs.pkl
python diff_cryodrgn_pkl.py output/pose.cs.pkl data/pose.cs.pkl

# Test write_starfile.py
cryodrgn write_starfile data/hand.5.mrcs output/ctf1.pkl -o output/test5.star

cryodrgn write_starfile data/hand.5.mrcs output/ctf1.pkl --ref-star data/FinalRefinement-OriginalParticles-PfCRT.star --keep-micrograph -o output/test6.star

echo All ok
