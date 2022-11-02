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

SCRIPT="$(pwd)/cryodrgn/commands/abinit_homo.py"
D128="$(pwd)/datasets/ribo_syn_128/projections.noise.wtrans.mrcs"
DREAL="$(pwd)/datasets/ribo_real_128/particles.128.phaseflip.mrcs"
DEFAULT="--t-extent 10 -n 5 --lr .0001 -b 8  --domain hartley --layers 3 --dim 256"

# N=fixedL_nkp4_bh2 run $SCRIPT $D128 $DEFAULT --nkeptposes 4 --base-healpy 2
# N=fixedL_nkp4_bh2_randinplane run $SCRIPT $D128 $DEFAULT --nkeptposes 4 --base-healpy 2
# N=fixedL_nkp24_bh2       run $SCRIPT $D128 $DEFAULT --nkeptposes 24 --base-healpy 2
# N=fixedL_nkp4_bh2_lmax31 run $SCRIPT $D128 $DEFAULT --nkeptposes 4  --base-healpy 2 --l-end 31
# N=fixedL_bh2_dim512      run $SCRIPT $D128 $DEFAULT --nkeptposes 4  --base-healpy 2 --dim 512
# N=fixedL_nkp4_bh2_n20    run $SCRIPT $D128 $DEFAULT --nkeptposes 4  --base-healpy 2 -n 20

# N=fixedL_nkp8_bh2_lmax48  run $SCRIPT $D128 $DEFAULT --nkeptposes 4 --base-healpy 2 --l-end 48
# N=fixedL_nkp4_bh2_lmax60  run $SCRIPT $D128 $DEFAULT --nkeptposes 4 --base-healpy 2 --l-end 60
# N=fixedL_nkp4_bh2_lmax60_lmin20  run $SCRIPT $D128 $DEFAULT --nkeptposes 4 --base-healpy 2 --l-end 48 --l-start 20 -b 4

# N=real1_end24             run $SCRIPT $DREAL --t-extent 20 --lr .0001 -b 8 --domain hartley --layers 3 --dim 256 --base-healpy 2 --nkeptposes 4 --niter 7 --l-end 24 -n 5
# N=real1_end48             run $SCRIPT $DREAL --t-extent 20 --lr .0001 -b 8 --domain hartley --layers 3 --dim 256 --base-healpy 2 --nkeptposes 4 --niter 7 --l-end 48 -n 5
# N=real1_end48_psf10       run $SCRIPT $DREAL --t-extent 20 --lr .0001 -b 8 --domain hartley --layers 3 --dim 256 --base-healpy 2 --nkeptposes 4 --niter 7 --l-end 48 -n 50 --ps-freq 10
# N=real1_end48_ramp5       run $SCRIPT $DREAL --t-extent 20 --lr .0001 -b 8 --domain hartley --layers 3 --dim 256 --base-healpy 2 --nkeptposes 4 --niter 7 --l-end 48 -n 10 --l-ramp-epochs 5
# N=real1_end48_ramp5_psf10 run $SCRIPT $DREAL --t-extent 20 --lr .0001 -b 8 --domain hartley --layers 3 --dim 256 --base-healpy 2 --nkeptposes 4 --niter 7 --l-end 48 -n 100 --ps-freq 10 --l-ramp-epochs 50
# N=real1_end24_ramp5       run $SCRIPT $DREAL --t-extent 20 --lr .0001 -b 8 --domain hartley --layers 3 --dim 256 --base-healpy 2 --nkeptposes 8 --niter 7 --l-end 24 -n 10 --l-ramp-epochs 5
# N=real1_end24_ramp5_psf10 run $SCRIPT $DREAL --t-extent 20 --lr .0001 -b 8 --domain hartley --layers 3 --dim 256 --base-healpy 2 --nkeptposes 8 --niter 7 --l-end 24 -n 100 --ps-freq 10 --l-ramp-epochs 50
# N=real1_end48_dim128        run $SCRIPT $DREAL --t-extent 20 --lr .0001 -b 8 --domain hartley --layers 3 --dim 128 --base-healpy 2 --nkeptposes 4 --niter 7 --l-end 48 -n 5
# N=real1_end48_dim128_psf10  run $SCRIPT $DREAL --t-extent 20 --lr .0001 -b 8 --domain hartley --layers 3 --dim 128 --base-healpy 2 --nkeptposes 4 --niter 7 --l-end 48 -n 50 --ps-freq 10
# N=real1_end48_rampP5       run $SCRIPT $DREAL --t-extent 20 --lr .0001 -b 8 --domain hartley --layers 3 --dim 256 --base-healpy 2 --nkeptposes 4 --niter 7 --l-end 48 -n 10 --l-ramp-epochs 5
# N=real1_end48_rampP5_psf10 run $SCRIPT $DREAL --t-extent 20 --lr .0001 -b 8 --domain hartley --layers 3 --dim 256 --base-healpy 2 --nkeptposes 4 --niter 7 --l-end 48 -n 100 --ps-freq 10 --l-ramp-epochs 50

SCRIPT="$(pwd)/cryodrgn/commands/abinit_het.py"

N=het1_zdim10_end48_rampP5_psf10_half_seed1 run $SCRIPT $DREAL --t-extent 20 --lr .0001 -b 8 --domain hartley --players 3 --pdim 256 --base-healpy 2 --nkeptposes 4 --niter 7 --l-end 48 -n 59 --ps-freq 10 --l-ramp-epochs 50 --zdim 10 --half-precision 1 --seed 1
N=het1_zdim10_end48_rampP5_psf10_half_seed2 run $SCRIPT $DREAL --t-extent 20 --lr .0001 -b 8 --domain hartley --players 3 --pdim 256 --base-healpy 2 --nkeptposes 4 --niter 7 --l-end 48 -n 59 --ps-freq 10 --l-ramp-epochs 50 --zdim 10 --half-precision 1 --seed 2
N=het1_zdim10_end48_rampP5_psf10_half_seed3 run $SCRIPT $DREAL --t-extent 20 --lr .0001 -b 8 --domain hartley --players 3 --pdim 256 --base-healpy 2 --nkeptposes 4 --niter 7 --l-end 48 -n 59 --ps-freq 10 --l-ramp-epochs 50 --zdim 10 --half-precision 1 --seed 3
