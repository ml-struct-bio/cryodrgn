set -e 
set -x

cryodrgn parse_ctf_star data/toy_projections.star -D 30 --Apix 1 -o test_ctf.pkl
python ../utils/write_starfile.py data/toy_projections.mrcs test_ctf.pkl -o output/test.star
python ../utils/write_starfile.py data/toy_projections.mrcs test_ctf.pkl -o output/test100.star --ind data/ind100.pkl 
