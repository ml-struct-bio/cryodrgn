set -e
set -x

# Test 1: Pase CTF from starfile with missing parameter
python ../utils/parse_ctf_star.py data/FinalRefinement-OriginalParticles-PfCRT.star -w .1 -o output/ctf1.pkl --png output/ctf1.png -N 5 -D 300 --Apix 1.035
diff output/ctf1.pkl data/ctf1.pkl

# Test 2: Parse poses from starfile
python ../utils/parse_pose_star.py data/FinalRefinement-OriginalParticles-PfCRT.star -o output/pose.star.pkl -D 300
diff output/pose.star.pkl data/pose.star.pkl

# Test 3: Parse poses from .cs file
python ../utils/parse_pose_csparc.py data/cryosparc_P12_J24_001_particles.cs -D 180 -o output/pose.cs.pkl
diff output/pose.cs.pkl data/pose.cs.pkl

echo All ok
