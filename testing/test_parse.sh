set -e
set -x

# Test 1: CTF parsing of starfile with missing parameter
python ../utils/parse_ctf_star.py data/FinalRefinement-OriginalParticles-PfCRT.star -w .1 -o output/ctf1.pkl --png output/ctf1.png -N 5 -D 300 --Apix 1.035
diff output/ctf1.pkl data/ctf1.pkl

echo All ok
