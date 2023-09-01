#!/bin/bash

set -e
set -x

cryodrgn train_vae data/sta_testing.star --datadir data/tilts/128 --encode-mode tilt --poses data/sta_pose.pkl --ctf data/sta_ctf.pkl --zdim 8 -o output/sta --tdim 256 --enc-dim 256 --dec-dim 256
