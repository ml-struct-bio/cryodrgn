python ../backproject.py data/toy_projections.mrcs data/toy_angles.pkl -o output/toy_recon.mrc
chimera output/toy_recon.mrc data/toymodel_small_nocenter.mrc 
python ../backproject_nn.py data/toy_projections.mrcs data/toy_angles.pkl -o output/toy_nn_recon -n 10
python ../backproject_hartley.py data/toy_projections.mrcs data/toy_angles.pkl -o output/toy_recon_hartley.mrc
python ../vae_rot.py  data/toy_projections.mrcs -o output/toy_recon_vae --lr .0001 
python ../vae_rot.py  data/aic_projections_small.mrc -o output/aic_projections --lr .0001 -n 10

python ../backproject_tilt.py /home/zhonge/research/cryoem/vae3d/00_data/hand/tilt_series/test/projections.mrcs /home/zhonge/research/cryoem/vae3d/00_data/hand/tilt_series/test/projections_tilt.mrcs /home/zhonge/research/cryoem/vae3d/00_data/hand/tilt_series/test/angles.pkl -o output/tilt/test.mrc --tilt -45
CUDA_VISIBLE_DEVICES=3 python ../vae_tilt.py /home/zhonge/research/cryoem/vae3d/00_data/hand/tilt_series/test/projections.mrcs /home/zhonge/research/cryoem/vae3d/00_data/hand/tilt_series/test/projections_tilt.mrcs -o output/tilt --tilt -45 --lr .0001
python ../project3d.py /home/zhonge/research/cryoem/vae3d/00_data/hand/hand_pointer.mrc -N 1000 -o output/project/projections.mrcs --out-rot output/project/rot.pkl --seed 0
python ../project3d.py /home/zhonge/research/cryoem/vae3d/00_data/hand/hand_pointer.mrc -N 1000 -o output/project/projections_tilt.mrcs --out-rot output/project/rot_tilt.pkl --seed 0 --tilt 45
python ../backproject.py output/project/projections.mrcs output/project/rot.pkl -o output/project/backproject.mrc --is-rot
python ../backproject_tilt.py output/project/projections.mrcs output/project/projections_tilt.mrcs output/project/rot.pkl --tilt 45 -o output/project/tilt_backproject.mrc --is-rot

python ../bnb_rot.py data/toy_projections.mrcs -o output/toy_recon_bnb --seed 0
python ../bnb_rot.py data/hand_projections.mrcs --load /red/zhonge/cryoem/vae3d/07_tilt/06_bnb/1/weights.pkl -o output/hand_bnb -n 11 --seed 82717
CUDA_VISIBLE_DEVICES=3 python ../vae_tilt.py data/tilt_series/projections.mrcs data/tilt_series/projections_tilt.mrcs --tilt 45 --lr .0005 -n 5 -o output/tilt > test1.log &

python ../backproject.py output/project/projections.mrcs output/project/rot.pkl --is-rot --trans output/project/trans.pkl -o output/project/backproject.mrc
python ../project3d.py /home/zhonge/research/cryoem/vae3d/00_data/hand/hand_pointer.mrc -N 1000 -o output/project/projections.mrcs --out-rot output/project/rot.pkl --seed 0 --out-trans output/project/trans.pkl
python ../vae_rot.py output/project/projections.mrcs --priors output/project/rot.pkl output/project/trans.pkl -o output/tmp
python ../vae_rot.py output/project/projections.mrcs --priors output/project/rot.pkl output/project/trans.pkl -o output/tmp
CUDA_VISIBLE_DEVICES=0 nohup python ../bnb_rot.py data/tilt_series/projections.mrcs -b 10 -n 2 -o output/bnb_trans > test1.log &
CUDA_VISIBLE_DEVICES=1 nohup python ../bnb_rot.py /home/zhonge/research/cryoem/vae3d/00_data/hand/tilt_series_wtrans/projections.mrcs -b 10 -n 2 -o output/bnb_trans2 > test2.log &
