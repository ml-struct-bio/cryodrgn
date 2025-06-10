#!/usr/bin/env sh
set -e
set -x

if [ $# -gt 2 ]; then
    echo "Usage: $0 [output_dir] [input_dir]"
    echo "  output_dir: Output directory (default: output)"
    echo "  input_dir: Input directory (default: data)"
    exit 1

elif [ $# -eq 0 ]; then
    OUTDIR=output
    INDIR=data
elif [ $# -eq 1 ]; then
    OUTDIR=$1
    INDIR=data
else
    INDIR=$1
    OUTDIR=$2
fi

particles=$INDIR/hand.mrcs
ctf=$INDIR/test_ctf.100.pkl
poses=$INDIR/hand_rot_trans.pkl


### Testing homogeneous reconstruction ###

# Fixed poses with cryoDRGN v3
cryodrgn train_nn $particles --poses $poses --ctf $ctf -o $OUTDIR/hand-recon_v3-fixed-hom --lr .0001 --seed 0

# Fixed poses with cryoDRGN-AI
cryodrgn setup $OUTDIR/hand-recon_v4-fixed-hom --particles $particles --poses $poses --ctf $ctf \
               --zdim 0 --model autodec --pose-estimation fixed
cryodrgn train $OUTDIR/hand-recon_v4-fixed-hom

# Ab-initio with cryoDRGN v3
cryodrgn abinit_homo $particles --ctf $ctf -o $OUTDIR/hand-recon_v3-abinit-hom --lr .0001 --seed 0
# Ab-initio with cryoDRGN-AI
cryodrgn setup $OUTDIR/hand-recon_v4-abinit-hom --particles $particles --ctf $ctf --zdim 0 --model autodec
cryodrgn train $OUTDIR/hand-recon_v4-abinit-hom

### Testing heterogeneous reconstruction ###

# Fixed poses with cryoDRGN v3
cryodrgn train_vae $particles --poses $poses --ctf $ctf -o $OUTDIR/hand-recon_v3-fixed-het \
                              --lr .0001 --seed 0 --zdim 10
cryodrgn analyze $OUTDIR/hand-recon_v3-fixed-het

# Fixed poses with cryoDRGN-AI
cryodrgn setup $OUTDIR/hand-recon_v4-fixed-het --particles $particles --poses $poses --ctf $ctf \
                                               --zdim 10 --model autodec --pose-estimation fixed
cryodrgn train $OUTDIR/hand-recon_v4-fixed-het

cryodrgn setup $OUTDIR/initial-conf-setup \
               --particles $particles --poses $poses --ctf $ctf \
               --zdim 10 --model autodec --pose-estimation fixed \
               --cfg load_z=$OUTDIR/hand-recon_v4-fixed-het/z.30.pkl lr_conf_table=0
cryodrgn train $OUTDIR/initial-conf-setup

exit 0

cryodrgn train $OUTDIR/hand-recon_v4-fixed-het_initial-conf \
               --from-outdir $OUTDIR/hand-recon_v4-fixed-het \
               --cfgs load_z=$OUTDIR/hand-recon_v4-fixed-het/z.10.pkl

cryodrgn train $OUTDIR/hand-recon_v4-fixed-het_initial-conf2 \
               --from-outdir $OUTDIR/hand-recon_v4-fixed-het \
               --cfgs load_z=$OUTDIR/hand-recon_v4-fixed-het/z.10.pkl lr_conf_table=0

# Ab-initio poses with cryoDRGN v3
cryodrgn abinit_het $particles --ctf $ctf -o $OUTDIR/hand-recon_v3-abinit-het --lr .0001 --seed 0 --zdim 10
cryodrgn analyze $OUTDIR/hand-recon_v3-abinit-het

# Ab-initio poses with cryoDRGN-AI
cryodrgn setup $OUTDIR/hand-recon_v4-abinit-het --particles $particles --ctf $ctf --zdim 10 --model autodec \
                                                --cfg t_ngrid=3 n_imgs_pose_search=100
cryodrgn train $OUTDIR/hand-recon_v4-abinit-het
