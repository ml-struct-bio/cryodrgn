#!/bin/bash
set -e

# https://stackoverflow.com/questions/59895
B=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

python $B/../cryodrgn/commands/abinit_het.py $B/data/hand.mrcs --zdim 8 -o $B/output/test --multigpu
python $B/../cryodrgn/commands/abinit_homo.py $B/data/hand.mrcs -o $B/output/test
