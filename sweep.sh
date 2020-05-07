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
# N=fixedL_nkp4_bh2_randinplane run $SCRIPT $D128 $DEFAULT --nkeptposes 4 --base-healpy 2
# N=fixedL_nkp24_bh2       run $SCRIPT $D128 $DEFAULT --nkeptposes 24 --base-healpy 2
N=fixedL_nkp4_bh2_lmax31 run $SCRIPT $D128 $DEFAULT --nkeptposes 4  --base-healpy 2 --l-end 31
# N=fixedL_bh2_dim512      run $SCRIPT $D128 $DEFAULT --nkeptposes 4  --base-healpy 2 --dim 512
# N=fixedL_nkp4_bh2_n20    run $SCRIPT $D128 $DEFAULT --nkeptposes 4  --base-healpy 2 -n 20

# N=fixedL_nkp8_bh2_lmax48  run $SCRIPT $D128 $DEFAULT --nkeptposes 4 --base-healpy 2 --l-end 48
# N=fixedL_nkp4_bh2_lmax60  run $SCRIPT $D128 $DEFAULT --nkeptposes 4 --base-healpy 2 --l-end 60
N=fixedL_nkp4_bh2_lmax60_lmin20  run $SCRIPT $D128 $DEFAULT --nkeptposes 4 --base-healpy 2 --l-end 48 --l-start 20 -b 4
