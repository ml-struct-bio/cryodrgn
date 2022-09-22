#!/bin/bash
set -e

B=$(dirname $0)

python $B/../cryodrgn/commands/abinit_het.py $B/data/hand.mrcs --zdim 8 -o $B/output/test
python $B/../cryodrgn/commands/abinit_homo.py $B/data/hand.mrcs -o $B/output/test
