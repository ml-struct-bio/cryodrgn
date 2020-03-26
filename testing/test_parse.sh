set -e
set -x

# Test 1: Pase CTF from starfile with missing parameter
python ../utils/parse_ctf_star.py data/FinalRefinement-OriginalParticles-PfCRT.star -w .1 -o output/ctf1.pkl --png output/ctf1.png -N 5 -D 300 --Apix 1.035
diff output/ctf1.pkl data/ctf1.pkl

# Test 2: Parse poses from starfile
python ../utils/parse_pose_star.py data/FinalRefinement-OriginalParticles-PfCRT.star -o output/parsed -D 300
diff output/parsed.rot.pkl data/parsed.rot.pkl
diff output/parsed.trans.pkl data/parsed.trans.pkl

# Test 3: Parse poses from .cs file
python ../utils/parse_pose_csparc.py data/cryosparc_P12_J24_001_particles.cs -D 180 -o output/parsed2
diff output/parsed2.rot.pkl data/parsed2.rot.pkl
diff output/parsed2.trans.pkl data/parsed2.trans.pkl

echo All ok
