#!/bin/bash
set -e

# https://stackoverflow.com/questions/59895
B=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

python $B/../cryodrgn/commands/abinit_het.py $B/data/hand.mrcs --zdim 8 -o $B/output/test --multigpu --domain hartley
python $B/../cryodrgn/commands/train_nn.py $B/data/hand.mrcs -o $B/output/test --domain hartley --uninvert-data --poses data/hand_rot.pkl --dim 256
python $B/../cryodrgn/commands/abinit_homo.py $B/data/hand.mrcs -o $B/output/test --load $B/output/test/weights.pkl --domain hartley -n 40 --uninvert-data
