#!/bin/bash

set -e

python setup.py develop

function run {
    O=/checkpoint/$USER/cryodrgn
    if [ -d $O/$N ]; then
        echo "Directory $O/$N exists; exiting"
        exit 1
    fi
    mkdir -p $O/$N
    git log -1 > $O/$N/GITLOG
    git diff >> $O/$N/GITLOG
    SETTINGS="-t 72:00:00 -J $N --partition dev \
              --output $O/$N/stdout.log --error $O/$N/stderr.log \
              --gres=gpu:1 --mem-per-gpu=64G --cpus-per-task 8 \
              --open-mode=append --chdir=$O/$N"

    OMP_NUM_THREADS=1 srun $SETTINGS -- python $@ -o $O/$N &
}

SCRIPT="$(pwd)/cryodrgn/commands/bnb_rot.py"
D128="$(pwd)/datasets/ribo_syn_128/projections.noise.wtrans.mrcs"
DEFAULT="--t-extent 10 -n 5 --lr .0001 -b 8  --domain hartley --layers 3 --dim 256"

# N=fixedL_nkp4_bh2 run $SCRIPT $D128 $DEFAULT --nkeptposes 4 --base-healpy 2
N=fixedL_nkp4_bh2_randinplane run $SCRIPT $D128 $DEFAULT --nkeptposes 4 --base-healpy 2
