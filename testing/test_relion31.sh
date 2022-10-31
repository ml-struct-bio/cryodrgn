set -e
set -x

# todo -- convert this starfile to 3.0 and compare outputs
cryodrgn downsample data/relion31.star -D 32 -o output/tmp.mrcs --relion31
cryodrgn downsample data/relion31.v2.star -D 32 -o output/tmp.mrcs --relion31
cryodrgn parse_pose_star data/relion31.star -D 256 -o output/pose.pkl --relion31 --Apix 1
cryodrgn parse_pose_star data/relion31.v2.star -D 256 -o output/pose.pkl --relion31 --Apix 1
cryodrgn parse_ctf_star data/relion31.star --Apix 1 -D 256 --relion31 --kv 300 -w .1 --ps 0 --cs 2.7 -o output/ctf.pkl
cryodrgn parse_ctf_star data/relion31.v2.star --Apix 1 -D 256 --relion31 --kv 300 -w .1 --ps 0 --cs 2.7 -o output/ctf.pkl
